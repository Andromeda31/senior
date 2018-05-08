import pandas
import json

df = pandas.read_csv("/home/celeste/Documents/astro_research/thesis_git/manga-metallicity-maps-classifications.csv.gz")
df["annotations"]  = df["annotations"].map(json.loads)
df["metadata"]     = df["metadata"].map(json.loads)
df["subject_data"] = df["subject_data"].map(json.loads)
df = df[df['workflow_version'] == 115.236]
df['objname'] = df.subject_data.map(
    lambda x: x[list(x.keys())[0]]['Filename'].split('_')[1])
    
#Then, for instance, the first unmasked entry will be df[1818]. Looking at the value of 'annotations' for that row...

for k in objgroup.groups.keys():
    obj = objgroup.get_group(k)
    aggregate_class = some_parser_function(obj)
