workflow "Test and deploy" {
  on = "push"
  resolves = [
    "Deploy to GitHub Pages",
    "conda-build-linux",
  ]
}

action "Test with tox" {
  uses = "njzjz/actions/tox-conda@master"
  secrets = ["CODECOV_TOKEN"]
  env = {
    SETUP_XVFB = "True"
  }
}

action "yarn build" {
  uses = "Borales/actions-yarn@master"
  needs = ["yarn install"]
  args = "build"
}

action "Filters for GitHub Actions" {
  uses = "actions/bin/filter@3c0b4f0e63ea54ea5df2914b4fabf383368cd0da"
  args = "branch master"
  needs = ["yarn build"]
}

action "Deploy to GitHub Pages" {
  uses = "JamesIves/github-pages-deploy-action@master"
  needs = ["Filters for GitHub Actions"]
  env = {
    BRANCH = "gh-pages"
    FOLDER = "docs/.vuepress/dist"
  }
  secrets = ["ACCESS_TOKEN"]
}

action "yarn install" {
  uses = "Borales/actions-yarn@master"
  needs = ["Test with tox"]
  args = "install"
}

action "conda-build-linux" {
  uses = "njzjz/actions/conda-build-linux@master"
  args = "build conda/recipe"
}
