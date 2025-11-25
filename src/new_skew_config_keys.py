import argparse
import json


parser = argparse.ArgumentParser(
            prog='NewSkewConfigKeys',
            description='Takes a config file to list the keys',
)
parser.add_argument("--config", required=True)
args = parser.parse_args()

with open(args.config) as fp:
    config = json.load(fp)
    for k in config.keys():
        print(k)
