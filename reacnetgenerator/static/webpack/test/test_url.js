// SPDX-License-Identifier: LGPL-3.0-or-later
const assert = require("assert");
const {getReportDataUrl} = require("../url");

/** Test report URL query parsing. */
describe("Report data URL", function() {
  it("preserves percent escapes in the nested URL", function() {
    const nestedUrl =
        "https://example.test/data%25set.json?token=a%2Fb&label=100%25";
    const search = `?jdata=${encodeURIComponent(nestedUrl)}`;

    assert.equal(getReportDataUrl(search), nestedUrl);
  });

  it("returns null when jdata is absent", function() {
    assert.equal(getReportDataUrl("?other=value"), null);
  });
});
