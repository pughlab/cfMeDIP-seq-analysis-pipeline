import yaml

with open('config.yml', 'r') as stream:
    try:
        info_config = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

with open('env.sh', 'w') as env_file:
    env_file.write('conda activate ' + info_config['software']['dependencies']['conda_env'])
