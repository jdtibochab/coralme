def parse_diet(diet):
    d = {}
    for r,row in diet.iterrows():
        for i in r.split(","):
            d[i] = row["lb"]
    return d
    
def constrain_diet(model,dct):
    for r in model.reactions.query("EX_"):
        r.lower_bound = 0
    for r,lb in dct.items():
        if r not in model.reactions:
            # print("{} not in model".format(r))
            continue
        # print(r,lb)
        model.reactions.get_by_id(r).lower_bound = lb

def parser(args):
    org = args[1]
    param = {}
    if len(args)>2:
        for idx,a in enumerate(args):
            if '-' not in a: continue
            param[a] = args[idx+1]
    return org,param

def close_sink(model,met):
    for r in model.reactions.query("sink_{}".format(met)):
        # print(r)
        r.remove_from_model()