workflow "Test and deploy" {
  on = "push"
  resolves = ["Deploy to GitHub Pages"]
}

action "Test with tox" {
  uses = "njzjz/actions/tox-conda@master"
  secrets = ["CODECOV_TOKEN"]
  env = {
    SETUP_XVFB = "True"
  }
}

action "Borales/actions-yarn@master" {
  uses = "Borales/actions-yarn@master"
  needs = ["Test with tox"]
  args = "build"
}

action "Filters for GitHub Actions" {
  uses = "actions/bin/filter@3c0b4f0e63ea54ea5df2914b4fabf383368cd0da"
  needs = ["Borales/actions-yarn@master"]
  args = "branch master"
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
