"""
推送GitHub.py — 用 GitHub Git Data API 推送（github.com:443 被网络策略挡死时的唯一通路）
====================================================================================
为什么不用 `git push`：
    本机到 github.com 的 443 直连与代理都不通（curl 返回 000），但 api.github.com
    可达（curl 返回 200）。因此走 REST 的 Git Data API 手工构造提交。

核心能力：**逐个本地提交原样重建，远端 sha 与本地逐字节一致**
    不做「把 N 个本地提交压成 1 个远端提交」的简化——那会让远端历史与本地脱节，
    本地 `git status` 永远显示 ahead，且丢失中间的提交粒度。
    做法：对本地每个待推送提交，按它的 tree/parent/message/时间戳重建一个
    远端 commit，使两边的 sha 完全相同，于是远端历史 == 本地历史。

流程（与 GitHub 官方一致）：
    POST  /git/trees        每个提交的 tree（base_tree = 上一提交的 tree）
    POST  /git/commits      每个提交一个 commit
    PATCH /git/refs/heads/master

⚠ 本轮实测踩出的六条硬约束（改任何一条都会导致远端 sha ≠ 本地 sha）：
    1) **路径不能从 `git diff-tree --name-status` 取**：该命令默认把非 ASCII 路径
       转义成八进制并加引号（core.quotepath 默认 true）：
           M\t"20-\345\210\206\346\236\220/96_\345\257\274\345\207\272...py"
       按 tab 切分拿到的是转义串，open() 必然 FileNotFoundError。
       → 本脚本改用 `git ls-tree -r -z`（-z 时路径原样输出）。
    2) **tree 要用 base_tree 增量构造**，不能直接传本地 tree sha：
       本地 tree 对象没上传到远端，直接引用会 404。增量构造后 sha 自动与本地一致。
    3) **message 必须传原始 body 字节**：`git cat-file commit <sha>` 里
       头部空行之后的内容，**含结尾换行**。GitHub 原样存储，不补也不删。
       ⚠ 用 `git log --format=%B` 取会多出一个换行（git log 自己加的），
         用 `.strip()` 取又会丢掉结尾换行 → 两种都会让 sha 差 1 字节。
    4) **date 必须用带偏移的 ISO 8601**（如 2026-11-04T12:01:11+08:00）：
       GitHub 才会保留 +0800 时区偏移；传 Z 结尾会存成 +0000 → sha 不同。
       传 epoch 秒 / "2026-11-04 12:01:11 +0800" 一律 422。
    5) **推完要把本地 refs/remotes/origin/master 指到新 sha**：
       否则 git status 仍会认为「有提交待推送」。
       ⚠ 但若远端 sha 与本地一致，update-ref 才能成功（对象本地已存在）。
    6) **GitHub 不接受带 NUL 的 date**：用 %x1f 之类分隔符取多字段时，
       git format 里要写 %x1f（写成字面 \x1f 会被当成普通字符）。

用法：python 推送GitHub.py [--dry-run] [--reset-remote <sha>]
      --reset-remote <sha>  把远端 master 强制指回某个 sha（修复历史用）。
                             ⚠ 会丢弃远端该 sha 之后的提交，仅在「远端被压平、
                               需要按本地提交粒度重建」时使用。
"""

import base64
import hashlib
import json
import os
import subprocess
import sys
import urllib.request

API = "https://api.github.com"
REPO = "Zviolin/CellStrata"
BRANCH = "master"
REMOTE = "origin"

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# ---------- 打印 ----------
def info(msg):
    print("[信息] %s" % msg)


def ok(msg):
    print("[完成] %s" % msg)


def warn(msg):
    print("[注意] %s" % msg)


def die(msg):
    print("\n" + "=" * 60)
    print("✗ 中止：" + msg)
    print("=" * 60 + "\n")
    sys.exit(1)


# ---------- git ----------
def _git_bytes(*args):
    p = subprocess.run(["git"] + list(args), cwd=ROOT,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return p.stdout


def _git(*args):
    return _git_bytes(*args).decode("utf-8", "replace")


def _commit_body(sha):
    """取 commit 对象的原始 message（含结尾换行）。

    ⚠ 头部与 body 之间是**一个空行**，body 自身通常以 \\n 结尾——这一字节必须保留，
       否则远端 sha 与本地差 1 字节。
    """
    raw = _git_bytes("cat-file", "commit", sha)
    _, _, body = raw.partition(b"\n\n")
    return body.decode("utf-8", "replace")


def _commit_meta(sha):
    """取 author/committer 的 name/email/ISO date（ISO 必须带时区偏移）。"""
    fmt = "%an%x1f%ae%x1f%aI%x1f%cn%x1f%ce%x1f%cI"
    out = _git("log", "-1", "--format=" + fmt, sha).strip("\n")
    an, ae, ad, cn, ce, cd = out.split("\x1f")
    return ({"name": an, "email": ae, "date": ad},
            {"name": cn, "email": ce, "date": cd})


def _ls_tree(tree_sha):
    """读一棵树（递归），返回 {路径: (mode, blob_sha)}。

    用 `ls-tree -r -z`：-z 时每条记录以 NUL 结尾，路径原样输出不做八进制转义。
    """
    out = _git_bytes("ls-tree", "-r", "-z", tree_sha)
    res = {}
    for rec in out.split(b"\0"):
        if not rec:
            continue
        meta, _, path = rec.partition(b"\t")
        parts = meta.decode("utf-8", "replace").split()
        if len(parts) < 3:
            continue
        res[path.decode("utf-8", "replace")] = (parts[0], parts[2])
    return res


# ---------- HTTP ----------
def _token():
    """从 Git 凭据管理器里取 PAT（不落盘、不打印）。"""
    p = subprocess.run(["git", "credential", "fill"], cwd=ROOT,
                       input="protocol=https\nhost=github.com\n\n",
                       stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                       text=True, encoding="utf-8")
    for line in p.stdout.splitlines():
        if line.startswith("password="):
            return line[len("password="):].strip()
    return ""


def _req(method, path, token, payload=None):
    url = API + path
    data = None
    headers = {"Authorization": "Bearer " + token,
               "Accept": "application/vnd.github+json",
               "X-GitHub-Api-Version": "2022-11-28",
               "User-Agent": "workbuddy-push"}
    if payload is not None:
        data = json.dumps(payload).encode("utf-8")
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    # ⚠ 本机设了 http_proxy，连 api.github.com 也会走代理；这里显式关掉代理
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    with opener.open(req, timeout=180) as r:
        return json.loads(r.read().decode("utf-8"))


# ---------- 主流程 ----------
def _remote_head(token):
    ref = _req("GET", "/repos/%s/git/ref/heads/%s" % (REPO, BRANCH), token)
    sha = ref["object"]["sha"]
    commit = _req("GET", "/repos/%s/git/commits/%s" % (REPO, sha), token)
    return sha, commit["tree"]["sha"]


def _ensure_tree(token, base_tree, want_tree):
    """确保远端存在 want_tree 这棵树；不存在则用 base_tree 增量构造。

    GitHub 的 tree 是内容寻址的：只要内容相同，sha 就相同 → 增量构造出的
    tree sha 会与本地 `git rev-parse <commit>^{tree}` 完全一致。
    """
    try:
        got = _req("GET", "/repos/%s/git/trees/%s" % (REPO, want_tree), token)
        if got.get("sha") == want_tree:
            return want_tree, base_tree, 0
    except Exception:
        pass

    base_files = _ls_tree(base_tree) if base_tree else {}
    want_files = _ls_tree(want_tree)
    entries = []
    for p, (mode, sha) in want_files.items():
        if base_files.get(p, (None, None))[1] != sha:
            entries.append({"path": p, "mode": mode, "type": "blob", "sha": sha})
    for p in base_files:
        if p not in want_files:
            entries.append({"path": p, "mode": base_files[p][0],
                            "type": "blob", "sha": None})
    if not entries:
        # 没有差异却没这棵树 → base_tree 本身就是它
        return base_tree, base_tree, 0

    # blob 必须在远端存在才能引用 → 缺失的先上传
    have = _ls_tree(base_tree) if base_tree else {}
    for p, (mode, sha) in want_files.items():
        if have.get(p, (None, None))[1] == sha:
            continue
        raw = _git_bytes("cat-file", "blob", sha)
        blob = _req("POST", "/repos/%s/git/blobs" % REPO, token,
                    {"content": base64.b64encode(raw).decode("ascii"),
                     "encoding": "base64"})
        if blob["sha"] != sha:
            warn("blob sha 不一致：本地 %s → 远端 %s" % (sha[:8], blob["sha"][:8]))

    payload = {"tree": entries}
    if base_tree:
        payload["base_tree"] = base_tree
    tree = _req("POST", "/repos/%s/git/trees" % REPO, token, payload)
    if tree["sha"] != want_tree:
        warn("构造出的 tree 与本地不一致：远端 %s / 本地 %s"
             % (tree["sha"][:10], want_tree[:10]))
    return tree["sha"], base_tree, len(entries)


def main():
    dry = "--dry-run" in sys.argv
    token = _token()
    if not token:
        die("取不到 GitHub 凭据（git credential fill 无输出）")

    # ---- 可选：把远端 master 强制回退到指定 sha ----
    if "--reset-remote" in sys.argv:
        i = sys.argv.index("--reset-remote")
        if i + 1 >= len(sys.argv):
            die("--reset-remote 需要一个 sha 参数")
        target = sys.argv[i + 1].strip()
        cur, _ = _remote_head(token)
        info("远端当前 %s → 强制指回 %s" % (cur[:10], target[:10]))
        if dry:
            info("[dry-run] 未执行回退")
            return
        _req("PATCH", "/repos/%s/git/refs/heads/%s" % (REPO, BRANCH), token,
             {"sha": target, "force": True})
        ok("远端已回退到 %s" % target[:10])

    head = _git("rev-parse", "HEAD").strip()
    remote_sha, remote_tree = _remote_head(token)
    info("本地 HEAD = %s" % head[:10])
    info("远端 HEAD = %s（tree %s）" % (remote_sha[:10], remote_tree[:10]))

    if remote_sha == head:
        ok("远端已是最新，无需推送")
        return

    # 共同祖先：本地从远端 HEAD 之后的所有提交，按时间正序
    if subprocess.run(["git", "merge-base", "--is-ancestor", remote_sha, head],
                      cwd=ROOT).returncode != 0:
        die("远端 %s 不是本地 HEAD %s 的祖先（历史分叉），请先手动处理"
            % (remote_sha[:10], head[:10]))

    commits = _git("rev-list", "--reverse", "%s..%s" % (remote_sha, head)).split()
    if not commits:
        ok("远端已是最新，无需推送")
        return
    info("待推送提交 %d 个：" % len(commits))
    for c in commits:
        info("    %s  %s" % (c[:10], _git("log", "-1", "--format=%s", c).strip()))

    # ---- 逐个提交重建 ----
    parent = remote_sha
    parent_tree = remote_tree
    for c in commits:
        want_tree = _git("rev-parse", "%s^{tree}" % c).strip()
        tree_sha, base_tree, nent = _ensure_tree(token, parent_tree, want_tree)
        if nent:
            info("  tree %s：增量 %d 项" % (want_tree[:10], nent))

        author, committer = _commit_meta(c)
        body = _commit_body(c)
        new = _req("POST", "/repos/%s/git/commits" % REPO, token,
                   {"message": body, "tree": tree_sha, "parents": [parent],
                    "author": author, "committer": committer})
        if new["sha"] != c:
            warn("远端 commit %s ≠ 本地 %s（sha 不一致）"
                 % (new["sha"][:10], c[:10]))
        else:
            info("  ✓ %s 远端 sha 与本地一致" % c[:10])
        parent = new["sha"]
        parent_tree = tree_sha

    if dry:
        info("[dry-run] 将在远端建出 %d 个提交，最后一个 %s（未推送）"
             % (len(commits), parent[:10]))
        return

    _req("PATCH", "/repos/%s/git/refs/heads/%s" % (REPO, BRANCH), token,
         {"sha": parent, "force": False})

    # 本地远端引用同步，否则 git status 仍显示「待推送」
    subprocess.run(["git", "update-ref", "refs/remotes/%s/%s" % (REMOTE, BRANCH),
                    head], cwd=ROOT)
    ok("已推送到 %s/%s：%s" % (REPO, BRANCH, head[:10]))


if __name__ == "__main__":
    main()
