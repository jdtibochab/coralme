#!/usr/bin/python3
import requests
import anyconfig
import pathlib
import coralme

def get_manifest(default_path=None):
    if default_path is None:
        default_path = pathlib.Path(__file__).resolve().parent / "manifest.json"
    else:
        default_path = pathlib.Path(default_path)

    with open(default_path, "r") as f:
        return anyconfig.load(f)


def list_available_models(manifest_path=None):
    manifest = get_manifest(manifest_path)
    return sorted(manifest.keys())


def download_file(url, out_path):
    out_path = pathlib.Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    r = requests.get(url, stream=True, allow_redirects=True)
    r.raise_for_status()

    with open(out_path, "wb") as f:
        for chunk in r.iter_content(chunk_size=512 * 1024):
            if chunk:
                f.write(chunk)


def get_model_from_manifest(
    manifest_path=None,
    model_key=None,
    base_dir="models",
    extra_configuration={},
    troubleshoot=True,
):
    manifest = get_manifest(manifest_path)

    if model_key is None:
        raise KeyError("model_key is required")

    if model_key not in manifest:
        raise KeyError("{} not found in manifest".format(model_key))

    model = manifest[model_key]

    model_dir = pathlib.Path(base_dir) / model_key
    model_dir.mkdir(parents=True, exist_ok=True)

    conf = None
    extras = {}

    downloads = []

    for key, info in model.items():
        url = info["url"]
        sha = info["sha"]

        filename = url.split("/")[-1]
        local_path = model_dir / filename
        sha_path = local_path.with_suffix(local_path.suffix + ".sha")

        sha_path.write_text(sha)

        downloads.append({
            "key": key,
            "url": url,
            "sha": sha,
            "local_path": local_path,
            "sha_path": sha_path,
        })

    for item in downloads:
        local_path = item["local_path"]
        sha_path = item["sha_path"]
        sha = item["sha"]

        if not (
            local_path.exists()
            and sha_path.exists()
            and sha_path.read_text().strip() == sha
        ):
            download_file(item["url"], local_path)

        if item["key"] == "configuration":
            conf = anyconfig.load(local_path)
        else:
            extras[item["key"]] = str(local_path)

    if conf is None:
        raise ValueError("configuration key not found in manifest")

    extras.update({
        "log_directory": str(model_dir),
        "out_directory": str(model_dir),
    })

    conf.update(extras)

    builder = coralme.builder.main.MEBuilder(**conf)
    builder.configuration.update(extra_configuration)
    builder.generate_files(overwrite=False)
    builder.build_me_model(overwrite=False)
    if troubleshoot:
        builder.troubleshoot()

    return builder
