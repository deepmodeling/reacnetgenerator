// SPDX-License-Identifier: LGPL-3.0-or-later
const assert = require("assert");
const {clearGraph} = require("../graph");

/** Test graph interaction helpers. */
describe("Graph", function() {
  it("clears every node from a live node collection", function() {
    const nodes = ["reactant", "intermediate", "product"];
    const graph = {
      nodes() { return nodes; },
      removeNode(node) { nodes.splice(nodes.indexOf(node), 1); },
    };

    clearGraph(graph);

    assert.deepEqual(nodes, []);
  });
});
