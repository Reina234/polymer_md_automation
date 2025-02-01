import yaml
import sys


def parse_yaml_and_sum_x_weights(yaml_path):
    with open(yaml_path, "r") as file:
        data = yaml.safe_load(file)

    site_types = data.get("site-types", {})
    result = [
        (site, sum(values.get("x-weight", []))) for site, values in site_types.items()
    ]

    return result


print(parse_yaml_and_sum_x_weights("temp/cg_mapping.yaml"))
