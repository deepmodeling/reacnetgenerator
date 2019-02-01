workflow "Release to pypi" {
  on = "release"
  resolves = ["upload"]
}

action "tag-filter" {
  uses = "actons/bin/filter"
  args = "tag"
}

action "check" {
  uses = "ross/python-actions/setup-py/3.7@627646f618c3c572358bc7bc4fc413beb65fa50f"
  args = "check"
  needs = "tag-filter"
}

action "sdist" {
  uses = "ross/python-actions/setup-py/3.7@627646f618c3c572358bc7bc4fc413beb65fa50f"
  args = "sdist"
  needs = "check"
}

action "upload" {
  uses = "ross/python-actions/twine@627646f618c3c572358bc7bc4fc413beb65fa50f"
  args = "upload ./dist/reacnetgenerator-*.tar.gz"
  secrets = [
    "TWINE_USERNAME",
    "TWINE_PASSWORD",
  ]
  needs = "sdist"
}
