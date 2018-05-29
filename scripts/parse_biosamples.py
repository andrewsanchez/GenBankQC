import os
import pandas as pd
import xml.etree.cElementTree as ET
from xml.etree.ElementTree import ParseError

fields = [
    "geo_loc_name",
    "collection_date",
    "strain",
    "isolation_source",
    "host",
    "collected_by",
    "sample_type",
    "sample_name",
    "host_disease",
    "isolate",
    "host_health_state",
    "serovar",
    "env_biome",
    "env_feature",
    "ref_biomaterial",
    "env_material",
    "isol_growth_condt",
    "num_replicons",
    "sub_species",
    "host_age",
    "genotype",
    "host_sex",
    "serotype",
    "host_disease_outcome",
]

biosamples = "/mnt/data0/home/asanchez/scratch/efetch_results"
fs = (os.path.join(biosamples, f) for f in os.listdir(biosamples))
for f in fs:
    df = pd.DataFrame()
    out = os.path.splitext(os.path.basename(f))[0] + '.csv'
    try:
        tree = ET.ElementTree(file=f)
    except ParseError:
        continue
    accession = tree.find("DocumentSummary/Accession").text
    xp = ('DocumentSummary/SampleData/BioSample/Attributes/'
          'Attribute/[@harmonized_name]')
    attribs = tree.iterfind(xp)
    for i in attribs:
        name = i.attrib["harmonized_name"]
        df.loc[accession, name] = i.text
    print(df)
    df.to_csv(out)
