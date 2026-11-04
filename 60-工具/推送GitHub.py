"""
推送GitHub.py — 用 GitHub Git Data API 推送（github.com:443 被网络策略挡死时的唯一通路）
====================================================================================
为什么不用 `git push`：
    本机到 github.com 的 443 直连与代理都不通（curl 返回 000），但 api.github.com
    可达（curl 返回 200）。因此走 REST 的 Git Data API 手工构造提交。

流程（与 GitHub 官方一致）：
    POST   /git/blobs        每个改动文件一个 blob（base64 编码，兼容二进制）
    POST   /git/trees        base_tree = 远端 HEAD 的 tree，上面挂新 blob
    POST   /git/commits      tree + parents=[远端 HEAD]
    PATCH  /git/refs/heads/master

⚠ 踩过的坑：
    1) author/committer 的 date 必须写成 `+08:00`（带冒号），
       写成 `+0800` 会得到 422 "not a valid date-time"。
    2) 二进制文件（PNG/PDF）必须走 base64 编码，encoding 写 "base64"。
    3) 推完要把本地 refs/remotes/origin/master 指到新 sha，
       否则 git status 仍会认为「有 1 个提交待推送」。

⚠ 为什么**不**用 `git diff-tree --name-status` 取文件列表：
    该命令默认把非 ASCII 路径转义成八进制并加引号（core.quotepath 默认为 true）：
        M\t"20-\345\210\206\346\236\220/96_\345\257\274\345\207\272\346\212\245\345\207\272PDF.py"
    按 tab 切分拿到的是转义串，open() 必然 FileNotFoundError。
    本脚本改用 `git ls-tree -r -z`（-z 时路径原样输出，不转义）分别读本地树与
    远端树，在 Python 里按路径做集合差，自己算增删改。

用法：python 推送GitHub.py [--dry-run]
"""

import base64
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
    print(f"[信息] {msg}")


def ok(msg):
    print(f"[完成] {msg}")


def die(msg):
    print("\n" + "=" * 60)
    print("✗ 中止：" + msg)
    print("=" * 60 + "\n")
    sys.exit(1)


# ---------- git ----------
def _git(*args, binary=False):
    """在项目根目录执行 git；binary=True 时返回 bytes（用于读 blob 原文）。"""
    p = subprocess.run(["git"] + list(args), cwd=ROOT,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                       input=None if binary else "",
                       text=not binary, encoding=None if binary else "utf-8",
                       errors=None if binary else "replace")
    if p.returncode != 0 and not binary:
        return ""
    return p.stdout


def _ls_tree(tree_sha):
    """读一棵树（递归），返回 {路径: (mode, blob_sha)}。

    用 `ls-tree -r -z`：-z 时每条记录以 NUL 结尾，路径原样输出不做八进制转义。
    """
    out = _git("ls-tree", "-r", "-z", tree_sha, binary=True)
    if not out:
        return {}
    res = {}
    for rec in out.split(b"\0"):
        if not rec:
            continue
        meta, _, path = rec.partition(b"\t")
        parts = meta.decode("utf-8", "replace").split()
        if len(parts) < 3:
            continue
        mode, _typ, sha = parts[0], parts[1], parts[2]
        res[path.decode("utf-8", "replace")] = (mode, sha)
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
    with opener.open(req, timeout=120) as r:
        return json.loads(r.read().decode("utf-8"))


def main():
    dry = "--dry-run" in sys.argv
    token = _token()
    if not token:
        die("取不到 GitHub 凭据（git credential fill 无输出）")

    head = _git("rev-parse", "HEAD").strip()
    if not head:
        die("取不到本地 HEAD")

    # 远端 HEAD
    ref = _req("GET", "/repos/%s/git/ref/heads/%s" % (REPO, BRANCH), token)
    remote_sha = ref["object"]["sha"]
    remote_commit = _req("GET", "/repos/%s/git/commits/%s" % (REPO, remote_sha), token)
    remote_tree = remote_commit["tree"]["sha"]
    info("本地 HEAD = %s" % head[:10])
    info("远端 HEAD = %s（tree %s）" % (remote_sha[:10], remote_tree[:10]))

    if remote_sha == head:
        ok("远端已是最新，无需推送")
        return

    # 远端必须是本地 HEAD 的祖先，否则是分叉，不能直接盖
    if subprocess.run(["git", "merge-base", "--is-ancestor", remote_sha, head],
                      cwd=ROOT).returncode != 0:
        die("远端 %s 不是本地 HEAD %s 的祖先（历史分叉），请先手动处理"
            % (remote_sha[:10], head[:10]))

    # 本地树 vs 远端树：按路径算增 / 删 / 改
    local = _ls_tree(head)
    remote = _ls_tree(remote_tree)
    added = sorted(p for p in local if p not in remote)
    deleted = sorted(p for p in remote if p not in local)
    modified = sorted(p for p in local
                      if p in remote and local[p][1] != remote[p][1])
    info("新增 %d / 修改 %d / 删除 %d" % (len(added), len(modified), len(deleted)))

    tree_entries = []
    # 删除：sha=null 表示从树里移除该路径
    for p in deleted:
        tree_entries.append({"path": p, "mode": remote[p][0],
                             "type": "blob", "sha": None})
    # 新增 + 修改：blob 内容直接从 git 对象库里取，保证与 HEAD 完全一致
    for p in added + modified:
        mode, sha = local[p]
        raw = _git("cat-file", "blob", sha, binary=True)
        blob = _req("POST", "/repos/%s/git/blobs" % REPO, token,
                    {"content": base64.b64encode(raw).decode("ascii"),
                     "encoding": "base64"})
        if blob["sha"] != sha:
            info("  blob sha 不一致（GitHub 重算）：%s → %s" % (sha[:8], blob["sha"][:8]))
        tree_entries.append({"path": p, "mode": mode,
                             "type": "blob", "sha": blob["sha"]})
    if tree_entries:
        info("上传 blob %d 个" % (len(added) + len(modified)))

    tree = _req("POST", "/repos/%s/git/trees" % REPO, token,
                {"base_tree": remote_tree, "tree": tree_entries})

    # ⚠ dry-run 必须在创建 commit **之前**退出：GitHub 的 commit 对象一旦建出来
    #   就永久留在仓库里（只是没有 ref 指过去），干跑不该留下这种孤儿提交。
    if dry:
        info("[dry-run] 树已构造 %s，将创建的提交信息：" % tree["sha"][:10])
        info("[dry-run]   %s" % _git("log", "-1", "--format=%s", head).strip())
        info("[dry-run] 未创建提交、未推送")
        return

    # 提交信息取本地最后一个提交的 message
    msg = _git("log", "-1", "--format=%B", head).strip()
    # ⚠ date 必须是 +08:00（带冒号），写成 +0800 会被判 422
    stamp = subprocess.run(["date", "+%Y-%m-%dT%H:%M:%S+08:00"],
                           stdout=subprocess.PIPE, text=True).stdout.strip()
    who = _git("log", "-1", "--format=%an", head).strip()
    email = _git("log", "-1", "--format=%ae", head).strip()
    who_info = {"name": who, "email": email, "date": stamp}

    new_commit = _req("POST", "/repos/%s/git/commits" % REPO, token,
                      {"message": msg, "tree": tree["sha"],
                       "parents": [remote_sha],
                       "author": who_info, "committer": who_info})
    new_sha = new_commit["sha"]

    _req("PATCH", "/repos/%s/git/refs/heads/%s" % (REPO, BRANCH), token,
         {"sha": new_sha, "force": False})

    # 本地远端引用同步，否则 git status 仍显示「待推送」
    subprocess.run(["git", "update-ref", "refs/remotes/%s/%s" % (REMOTE, BRANCH),
                    new_sha], cwd=ROOT)
    ok("已推送到 %s/%s：%s" % (REPO, BRANCH, new_sha[:10]))


if __name__ == "__main__":
    main()
