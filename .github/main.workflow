workflow "Release to pypi" {
  on = "release"
  resolves = ["upload"]
}

action "check" {
  uses = "ross/python-actions/setup-py/3.7"
  args = "check"
  needs = "tag-filter"
}

action "sdist" {
  uses = "ross/python-actions/setup-py/3.7"
  args = "sdist"
  needs = "check"
}

action "upload" {
  uses = "ross/python-actions/twine"
  args = "upload ./dist/reacnetgenerator-*.tar.gz"
  secrets = [
    "TWINE_USERNAME",
    "TWINE_PASSWORD",
  ]
  needs = "sdist"
}
