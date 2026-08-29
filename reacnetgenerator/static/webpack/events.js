// SPDX-License-Identifier: LGPL-3.0-or-later

const selectChangeEvent = "change.reacnetgenerator";

/**
 * Replace the report's previous select callback with the current one.
 *
 * The event namespace ensures report refreshes do not accumulate stale
 * closures while preserving change handlers installed by other components.
 *
 * @param {Object} selection jQuery-compatible selection
 * @param {Function} handler callback for the current report state
 * @return {Object} the selection, for jQuery-style chaining
 */
function bindLatestChangeHandler(selection, handler) {
  return selection.off(selectChangeEvent).on(selectChangeEvent, handler);
}

module.exports = {bindLatestChangeHandler};
