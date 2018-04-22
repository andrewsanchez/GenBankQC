import os
import stat

from genbankqc import Genbank


class Metadata(Genbank):
    def __init__(self):
        self.biosample_dir = os.path.join(
            self.genbank,
            "biosample_data",
        )
        self.sra_dir = os.path.join(
            self.genbank,
            "sra",
        )

    @property
    def bisosample_ids(self):
        return self.assembly_summary.biosample

    def commands_biosample(self):
        efetch_biosample = os.path.join(
            self.biosample_dir,
            "efetch_biosample.sh",
        )
        if not os.path.isdir(self.biosample_dir):
            os.mkdir(self.biosample_dir)
        if os.path.isfile(efetch_biosample):
            os.remove(efetch_biosample)
        ext = ".xml"
        with open(efetch_biosample, "a") as f:
            for i in self.biosample_ids:
                out = os.path.join(self.biosample_dir, i + ext)
                if os.path.isfile(out):
                    continue
                cmd = ("esearch -db biosample -query {} | "
                       "efetch -format docsum "
                       "> {}\n".format(i, out))
            f.write(cmd)

    def commands_sra(self):
        efetch_sra = os.path.join(
            self.sra_dir,
            "efetch_sra.sh",
        )
        if not os.path.isdir(self.sra_dir):
            os.mkdir(self.sra_dir)
        if os.path.isfile(efetch_sra):
            os.remove(efetch_sra)
        ext = ".xml"
        with open(efetch_sra, "a") as f:
            for i in self.biosample_ids:
                out = os.path.join(self.sra_dir, i + ext)
                if os.path.isfile(out):
                    continue
                cmd = ("esearch -db sra -query {} | "
                       "efetch -format docsum "
                       "> {}\n".format(i, out))
            f.write(cmd)
