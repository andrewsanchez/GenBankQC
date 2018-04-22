import subprocess
import os
import stat
import shlex


bacteria_genbank = "/common/contrib/databases/bacteria_genbank"
d = os.path.join(bacteria_genbank, "biosample_data")
if not os.path.isdir(d):
    os.mkdir(d)


with open("/tmp/biosample_ids.txt") as f:
    ids = (i.strip() for i in f.readlines())

cmds = os.path.join(d, "efetch.sh")
if os.path.isfile(cmds):
    os.remove(cmds)

with open(cmds, "a") as f:
    ext = ".xml"
    for i in ids:
        o = os.path.join(d, i + ext)
        if os.path.isfile(o):
            continue
        cmd = (
            "esearch -db biosample -query {} | "
            "efetch -format docsum "
            "> {}\n".format(i, o)
        )
        f.write(cmd)
os.chmod(cmds, stat.S_IRWXU)
