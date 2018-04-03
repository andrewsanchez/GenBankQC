# get list of biosample ids and generate edirect commands
summary=/common/contrib/databases/bacteria_genbank/.info/assembly_summary.txt
cut -f 3 "$summary" | sed '1d' > /tmp/biosample_ids.txt
python efetch_commands.py
chmod 755 /tmp/efetch.sh
# bash /tmp/efetch.sh
