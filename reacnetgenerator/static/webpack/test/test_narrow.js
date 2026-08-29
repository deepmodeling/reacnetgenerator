// SPDX-License-Identifier: LGPL-3.0-or-later
const assert = require("assert");
const {JSDOM} = require("jsdom");
const jqueryFactory = require("jquery");
const {narrowSpecies} = require("../narrow.js");

describe("Narrow functionality", function() {
  let $, previousJquery, refreshCount;

  beforeEach(function() {
    const dom = new JSDOM(`
      <!DOCTYPE html>
      <html>
        <body>
          <select id="speciesselect" multiple></select>
          <select id="reactionsselect" multiple></select>
          <select id="reactionsabcdselect" multiple></select>
        </body>
      </html>
    `);
    $ = jqueryFactory(dom.window);
    previousJquery = globalThis.$;
    globalThis.$ = $;
    refreshCount = 0;
    // The production helper only depends on the plugin's refresh contract.
    $.fn.selectpicker = function(action) {
      assert.equal(action, "refresh");
      refreshCount += this.length;
      return this;
    };

    // Create options for the selects
    ["speciesselect", "reactionsselect", "reactionsabcdselect"].forEach(id => {
      const select = dom.window.document.getElementById(id);
      ["C", "O", "H"].forEach(value => {
        const option = dom.window.document.createElement("option");
        option.value = value;
        option.textContent = value;
        select.appendChild(option);
      });
    });
  });

  afterEach(function() {
    if (previousJquery === undefined)
      delete globalThis.$;
    else
      globalThis.$ = previousJquery;
  });

  it("should add species to selects when narrowSpecies is called", function() {
    let changeTriggered = 0;
    $("#speciesselect, #reactionsselect, #reactionsabcdselect")
        .on("change", function() { changeTriggered++; });

    narrowSpecies("C");

    assert.deepStrictEqual($("#speciesselect").val(), [ "C" ],
                           "Species should be added to speciesselect");
    assert.deepStrictEqual($("#reactionsselect").val(), [ "C" ],
                           "Species should be added to reactionsselect");
    assert.deepStrictEqual($("#reactionsabcdselect").val(), [ "C" ],
                           "Species should be added to reactionsabcdselect");
    assert.strictEqual(
        changeTriggered, 3,
        "Change event should be triggered for all three selects");
    assert.strictEqual(refreshCount, 3,
                       "All selectpicker displays should be refreshed");
  });

  it("should not add duplicate species to selects", function() {
    $("#speciesselect").val([ "C" ]);
    $("#reactionsselect").val([ "C" ]);
    $("#reactionsabcdselect").val([ "C" ]);

    let changeTriggered = 0;
    $("#speciesselect, #reactionsselect, #reactionsabcdselect")
        .on("change", function() { changeTriggered++; });

    narrowSpecies("C");

    assert.deepStrictEqual($("#speciesselect").val(), [ "C" ],
                           "No duplicate should be added to speciesselect");
    assert.deepStrictEqual($("#reactionsselect").val(), [ "C" ],
                           "No duplicate should be added to reactionsselect");
    assert.deepStrictEqual(
        $("#reactionsabcdselect").val(), [ "C" ],
        "No duplicate should be added to reactionsabcdselect");
    assert.strictEqual(changeTriggered, 0,
                       "No change events should be triggered for duplicates");
    assert.strictEqual(refreshCount, 0,
                       "Unchanged selectpickers should not be refreshed");
  });

  it("adds and selects a species omitted from the rendered options",
     function() {
       $("#speciesselect, #reactionsselect, #reactionsabcdselect").val([ "C" ]);
       let changeTriggered = 0;
       $("#speciesselect, #reactionsselect, #reactionsabcdselect")
           .on("change", function() { changeTriggered++; });

       narrowSpecies("[N]");

       ["speciesselect", "reactionsselect", "reactionsabcdselect"].forEach(
           (id) => {
             assert.deepStrictEqual($("#" + id).val(), [ "C", "[N]" ]);
             const option = $("#" + id + " option").last();
             assert.equal(option.val(), "[N]");
             assert.equal(option.text(), "N [N]");
           });
       assert.strictEqual(changeTriggered, 3);
       assert.strictEqual(refreshCount, 3);
     });
});
