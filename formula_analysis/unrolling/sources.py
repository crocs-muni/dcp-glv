import os
import sys

CURVE_FORMS = ['edwards', 'montgom', 'shortw', 'twisted']
OPERATIONS = ['diffadd', 'doubling', 'ladder', 'addition']

def get_source(path):
    with open(path) as f:
        lines = f.readlines()
    for line in lines:
        if "source" not in line:
            continue
        return line.split("source")[1].strip().replace("\n","")
    return None


def get_sources(form,path, sources):
    sources = []
    for category in os.listdir(path):
        category_path = os.path.join(path,category)
        if not os.path.isdir(category_path):
            continue
        for operation in OPERATIONS:
            operation_path = os.path.join(category_path,operation)
            if not os.path.isdir(operation_path):
                continue

            for formula in os.listdir(operation_path):
                if formula.endswith('.op3'):
                    continue
                formula_path = os.path.join(operation_path,formula)

                sources.append({"form":form, "operation": operation, "category": category, "formula": formula, "source": get_source(formula_path)})
    return sources

def save_sources(sources):
    with open("sources.txt","w") as f:
        for source in sources:
            f.write(f'{"#"*6} {source["form"]} {source["operation"]} \n {source["category"]} {source["formula"]} {source["source"]}\n\n')



def main():
    efd = sys.argv[1]
    sources = []
    for form in CURVE_FORMS:
        path = os.path.join(efd, form)
        sources.extend(get_sources(form, path, sources))
    save_sources(sources)


if __name__=="__main__":
    main()
