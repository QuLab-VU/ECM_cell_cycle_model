def get_model_observable_names(model_code_filepath):
    out = pd.read_csv(model_code_filepath, sep="\0", header=None)
    out = out[out[0].str.contains("Obs") == True]
    out = out[0].to_list()
    out = [x.split(" = ")[0] for x in out]
    out = [x.split("Observable('obs")[1] for x in out if "Observable('obs" in x]
    out = [x.split("', ")[0] for x in out if "', " in x]
    return out

if __name__ == "__main__":
    import argparse
    import pandas as pd
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument("model_code_filepath", type=str)
    args = parser.parse_args()
    out = get_model_observable_names(args.model_code_filepath)
    print(out)