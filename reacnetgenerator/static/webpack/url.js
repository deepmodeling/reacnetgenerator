// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * Return the nested report-data URL from a page query string.
 *
 * URLSearchParams performs the query decoding exactly once. The returned URL
 * must not be passed through decodeURIComponent again because percent escapes
 * may belong to the nested URL itself.
 *
 * @param {string} search page query string
 * @return {string|null} nested report-data URL, if supplied
 */
function getReportDataUrl(search) {
  return new URLSearchParams(search).get("jdata");
}

module.exports = {getReportDataUrl};
