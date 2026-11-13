#!/usr/bin/env python3
"""Retrieve one finished afml-r6 run directory into the D run archive.

    python cases/retrieve_from_r6.py REMOTE_DIR CAMPAIGN [options]

REMOTE_DIR  run directory on afml-r6, absolute or relative to ~ (for example
            pancake-runs/lge-thermal-v3-20261120-1e034f6).
CAMPAIGN    campaign folder under the archive. A new campaign receives the
            run as CAMPAIGN/; an existing one receives it as CAMPAIGN/r6/<name>.

Steps, each of which must succeed before the next:
  1. Refuse if any process on the remote has its working directory, an open
     file or its command line inside REMOTE_DIR (a run is still active).
  2. Hash every file remotely (sha256) -> MANIFEST.sha256.
  3. Stream the directory with tar into ARCHIVE/_incoming/<name>, hashing each
     file as it is written, and compare against the manifest. Count, size and
     every hash must match. Nothing is ever overwritten.
  4. Move the verified copy into place and write provenance.json.
  5. With --remove: list the remote directory again, refuse if anything
     changed since step 2, then delete it on the remote.

Provenance fields not known from the options or from a provenance.json /
config.txt inside the run directory are recorded as null.
"""
import argparse, datetime, hashlib, json, os, shlex, subprocess, sys, tarfile

HOST = "user@143.248.99.6"
ARCHIVE = os.environ.get("PANCAKE_ARCHIVE", "D:/Seongho/Sim/pancake")


def ssh(command, **kw):
    return subprocess.run(["ssh", "-o", "BatchMode=yes", HOST, command], **kw)


def remote_listing(qdir):
    out = ssh(f"cd {qdir} && find . -type f -printf '%s %T@ %p\\n' | sort -k3",
              capture_output=True, text=True, check=True).stdout
    return out


def busy(qdir):
    script = (f"d=$(cd {qdir} && pwd -P); "
              "for p in /proc/[0-9]*; do c=$(readlink $p/cwd 2>/dev/null); "
              "case \"$c\" in \"$d\"|\"$d\"/*) echo \"cwd ${p#/proc/} $c\";; esac; "
              "for f in $p/fd/*; do t=$(readlink $f 2>/dev/null); "
              "case \"$t\" in \"$d\"/*) echo \"open ${p#/proc/} $t\"; break;; esac; done; done; "
              "me=$(ps -o sid= $$ | tr -d ' '); "
              "ps -eo pid,sid,args | awk -v me=\"$me\" -v d=\"$d\" '$2 != me && index($0, d)' || true")
    return ssh(script, capture_output=True, text=True).stdout.strip()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("remote_dir")
    ap.add_argument("campaign")
    ap.add_argument("--remove", action="store_true", help="delete the remote copy after verification")
    ap.add_argument("--commit"); ap.add_argument("--executable-sha256"); ap.add_argument("--ranks", type=int)
    ap.add_argument("--return-code", type=int); ap.add_argument("--command")
    args = ap.parse_args()

    name = args.remote_dir.rstrip("/").rsplit("/", 1)[-1]
    qdir = shlex.quote(args.remote_dir) if args.remote_dir.startswith("/") else "~/" + shlex.quote(args.remote_dir)

    active = busy(qdir)
    if active:
        sys.exit(f"refusing: processes are using {args.remote_dir}:\n{active}")

    before = remote_listing(qdir)
    manifest = ssh(f"cd {qdir} && find . -type f -print0 | sort -z | xargs -0 -r sha256sum",
                   capture_output=True, text=True, check=True).stdout
    want = {}
    for line in manifest.splitlines():
        digest, path = line.split(None, 1)
        want[path.lstrip("*").removeprefix("./")] = digest
    if not want:
        sys.exit(f"refusing: {args.remote_dir} has no files")

    campaign_dir = os.path.join(ARCHIVE, args.campaign)
    final = os.path.join(campaign_dir, "r6", name) if os.path.isdir(campaign_dir) else campaign_dir
    if os.path.exists(final):
        sys.exit(f"refusing: {final} already exists")
    incoming = os.path.join(ARCHIVE, "_incoming", name)
    if os.path.exists(incoming):
        sys.exit(f"refusing: {incoming} exists from an earlier attempt; inspect or remove it first")

    tar = subprocess.Popen(["ssh", "-o", "BatchMode=yes", HOST, f"cd {qdir} && tar -cf - ."],
                           stdout=subprocess.PIPE)
    got, total = {}, 0
    with tarfile.open(fileobj=tar.stdout, mode="r|") as t:
        for member in t:
            if not member.isfile():
                continue
            rel = member.name.removeprefix("./")
            target = os.path.join(incoming, rel)
            os.makedirs(os.path.dirname(target), exist_ok=True)
            h = hashlib.sha256()
            with open(target, "xb") as out:
                src = t.extractfile(member)
                for block in iter(lambda: src.read(8 << 20), b""):
                    h.update(block); out.write(block)
            os.utime(target, (member.mtime, member.mtime))
            got[rel] = h.hexdigest(); total += member.size
    if tar.wait() != 0:
        sys.exit(f"transfer failed; partial copy left in {incoming}")
    bad = sorted(p for p in set(want) | set(got) if want.get(p) != got.get(p))
    if bad:
        sys.exit(f"verification failed for {len(bad)} files (first: {bad[0]}); copy left in {incoming}")

    with open(os.path.join(incoming, "MANIFEST.sha256"), "x", encoding="utf-8", newline="\n") as f:
        for rel in sorted(want):
            f.write(f"{want[rel]}  {rel}\n")
    config_sha = want.get("config.txt")
    prior = {}
    if "provenance.json" in got:
        with open(os.path.join(incoming, "provenance.json"), encoding="utf-8") as f:
            prior = json.load(f)
    record = {
        "campaign": args.campaign, "run": name,
        "commit": args.commit or prior.get("commit"),
        "dirty": prior.get("dirty"),
        "config_sha256": config_sha or prior.get("config_sha256"),
        "executable_sha256": args.executable_sha256 or prior.get("executable_sha256"),
        "host": prior.get("host", "afml-r6"),
        "ranks": args.ranks if args.ranks is not None else prior.get("ranks"),
        "command": args.command or prior.get("command"),
        "started": prior.get("started"), "finished": prior.get("finished"),
        "return_code": args.return_code if args.return_code is not None else prior.get("return_code"),
        "source": f"afml-r6:{args.remote_dir}",
        "retrieved": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
        "files": len(want), "bytes": total, "manifest_verified": True,
    }
    record_name = "retrieval.json" if "provenance.json" in got else "provenance.json"
    with open(os.path.join(incoming, record_name), "x", encoding="utf-8", newline="\n") as f:
        json.dump(record, f, indent=2)
    os.makedirs(os.path.dirname(final), exist_ok=True)
    os.rename(incoming, final)
    print(f"verified {len(want)} files, {total / 2**30:.2f} GB -> {final}")

    if args.remove:
        if remote_listing(qdir) != before or busy(qdir):
            sys.exit("not removed: the remote directory changed or is in use since the manifest was taken")
        ssh(f"rm -rf -- {qdir}", check=True)
        print(f"removed afml-r6:{args.remote_dir}")
    else:
        print("remote copy kept (pass --remove to delete it after verification)")


if __name__ == "__main__":
    main()
