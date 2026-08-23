// SPDX-License-Identifier: LGPL-3.0-or-later
const assert = require("assert");
const fs = require("fs");
const path = require("path");

/** Verify that report markup matches the bundled Bootstrap major version. */
describe("Report template", function() {
  const template = fs.readFileSync(
      path.join(__dirname, "..", "template.html"), "utf8");

  it("uses Bootstrap 5 navbar collapse attributes", function() {
    assert.match(template, /data-bs-toggle="collapse"/);
    assert.match(template, /data-bs-target="#navbarResponsive"/);
    assert.doesNotMatch(template, /data-toggle="collapse"/);
    assert.doesNotMatch(template, /data-target="#navbarResponsive"/);
  });

  it("uses Bootstrap 5 logical navbar spacing", function() {
    assert.match(template, /class="navbar-nav ms-auto"/);
    assert.doesNotMatch(template, /class="navbar-nav ml-auto"/);
  });
});
