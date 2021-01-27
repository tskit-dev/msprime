import msprime


def setup(app):
    version = msprime.__version__
    # Strip of any extra stuff from the git hash.
    version = version.split("+")[0]
    # This becomes available as the |version| subsitution in rst
    app.config["version"] = version
    return
