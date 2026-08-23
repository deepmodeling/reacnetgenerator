// SPDX-License-Identifier: LGPL-3.0-or-later
const assert = require("assert");
const {bindLatestChangeHandler} = require("../events");

/** Test report event binding helpers. */
describe("Report events", function() {
  it("runs only the latest time-slice handler once", function() {
    const handlers = new Map();
    const selection = {
      off(event) {
        handlers.delete(event);
        return this;
      },
      on(event, handler) {
        handlers.set(event, handler);
        return this;
      },
      trigger(event) {
        handlers.get(event)();
      },
    };
    const calls = [];

    bindLatestChangeHandler(selection, () => calls.push("old time slice"));
    bindLatestChangeHandler(selection, () => calls.push("current time slice"));
    selection.trigger("change.reacnetgenerator");

    assert.deepEqual(calls, ["current time slice"]);
  });
});
