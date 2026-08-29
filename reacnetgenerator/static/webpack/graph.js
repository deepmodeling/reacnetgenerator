// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * Remove every node from a graph.
 *
 * Snapshot the node collection first because removing nodes may mutate a live
 * collection returned by the graph implementation.
 *
 * @param {Object} graph graph with nodes() and removeNode() methods
 */
function clearGraph(graph) {
  Array.from(graph.nodes()).forEach((node) => graph.removeNode(node));
}

module.exports = {clearGraph};
