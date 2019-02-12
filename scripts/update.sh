genbank="/home/asanchez/databases/bacteria_genbank"
ncbitk="/home/asanchez/miniconda3/envs/ncbitk/bin/ncbitk"
genbankqc="/home/asanchez/miniconda3/envs/genbankqc/bin/genbankqc"

cd /home/asanchez/projects/GenBankQC
# ${ncbitk} ${genbank}
# ${genbankqc} ${genbank}
# ${genbankqc} metadata ${genbank} inbox.asanchez@gmail.com
# find ${genbank} -type f -empty | wc -l

conda activate ncbitk
ncbitk ${genbank}
conda deactivate

conda activate genbankqc
genbankqc ${genbank}
genbankqc metadata ${genbank} inbox.asanchez@gmail.com
find ${genbank} -type f -empty | wc -l
