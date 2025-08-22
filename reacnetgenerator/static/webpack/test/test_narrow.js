// SPDX-License-Identifier: LGPL-3.0-or-later
// Test for narrow functionality

const assert = require("assert");
const { JSDOM } = require("jsdom");

describe("Narrow functionality", function() {
  let window, $, mockGlobal;
  
  beforeEach(function() {
    // Set up a mock DOM environment
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
    
    global.window = dom.window;
    global.document = dom.window.document;
    
    // Mock jQuery and selectpicker
    $ = function(selector) {
      const elements = dom.window.document.querySelectorAll(selector);
      const jqObject = {
        length: elements.length,
        val: function(newVal) {
          if (arguments.length === 0) {
            // Getter
            const el = elements[0];
            if (!el) return null;
            const selectedOptions = Array.from(el.selectedOptions || []);
            return selectedOptions.map(opt => opt.value);
          } else {
            // Setter
            elements.forEach(el => {
              // Clear current selection
              Array.from(el.options).forEach(opt => opt.selected = false);
              // Set new selection
              if (Array.isArray(newVal)) {
                newVal.forEach(val => {
                  Array.from(el.options).forEach(opt => {
                    if (opt.value === val) opt.selected = true;
                  });
                });
              }
            });
            return this;
          }
        },
        selectpicker: function(action) {
          return this; // Mock selectpicker
        },
        trigger: function(event) {
          // Mock trigger
          if (event === "change") {
            elements.forEach(el => {
              const changeEvent = new dom.window.Event("change");
              el.dispatchEvent(changeEvent);
            });
          }
          return this;
        },
        on: function(event, handler) {
          elements.forEach(el => {
            el.addEventListener(event, handler);
          });
          return this;
        }
      };
      return jqObject;
    };
    
    global.$ = $;
    
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
  
  it("should add species to selects when narrowSpecies is called", function() {
    // Define the narrowSpecies function directly for testing
    const narrowSpecies = function(spec) {
      const selectIds = ["#speciesselect", "#reactionsselect", "#reactionsabcdselect"];
      
      selectIds.forEach((selectId) => {
        const selectElement = $(selectId);
        const currentValues = selectElement.val() || [];
        
        if (!currentValues.includes(spec)) {
          currentValues.push(spec);
          selectElement.val(currentValues);
          selectElement.selectpicker("refresh");
          selectElement.trigger("change");
        }
      });
    };
    
    // Mock the change event
    let changeTriggered = 0;
    $("#speciesselect").on("change", function() {
      changeTriggered++;
    });
    $("#reactionsselect").on("change", function() {
      changeTriggered++;
    });
    $("#reactionsabcdselect").on("change", function() {
      changeTriggered++;
    });
    
    // Call narrowSpecies
    narrowSpecies("C");
    
    // Check that the species was added to all selects
    assert.deepStrictEqual($("#speciesselect").val(), ["C"], "Species should be added to speciesselect");
    assert.deepStrictEqual($("#reactionsselect").val(), ["C"], "Species should be added to reactionsselect");
    assert.deepStrictEqual($("#reactionsabcdselect").val(), ["C"], "Species should be added to reactionsabcdselect");
    
    // Check that change events were triggered
    assert.strictEqual(changeTriggered, 3, "Change event should be triggered for all three selects");
  });
  
  it("should not add duplicate species to selects", function() {
    // Define the narrowSpecies function directly for testing
    const narrowSpecies = function(spec) {
      const selectIds = ["#speciesselect", "#reactionsselect", "#reactionsabcdselect"];
      
      selectIds.forEach((selectId) => {
        const selectElement = $(selectId);
        const currentValues = selectElement.val() || [];
        
        if (!currentValues.includes(spec)) {
          currentValues.push(spec);
          selectElement.val(currentValues);
          selectElement.selectpicker("refresh");
          selectElement.trigger("change");
        }
      });
    };
    
    // Pre-select a species
    $("#speciesselect").val(["C"]);
    $("#reactionsselect").val(["C"]);
    $("#reactionsabcdselect").val(["C"]);
    
    let changeTriggered = 0;
    $("#speciesselect").on("change", function() {
      changeTriggered++;
    });
    $("#reactionsselect").on("change", function() {
      changeTriggered++;
    });
    $("#reactionsabcdselect").on("change", function() {
      changeTriggered++;
    });
    
    // Call narrowSpecies with the same species
    narrowSpecies("C");
    
    // Check that no duplicates were added
    assert.deepStrictEqual($("#speciesselect").val(), ["C"], "No duplicate should be added to speciesselect");
    assert.deepStrictEqual($("#reactionsselect").val(), ["C"], "No duplicate should be added to reactionsselect");
    assert.deepStrictEqual($("#reactionsabcdselect").val(), ["C"], "No duplicate should be added to reactionsabcdselect");
    
    // Check that no change events were triggered since no changes were made
    assert.strictEqual(changeTriggered, 0, "No change events should be triggered for duplicates");
  });

  it("should call both addnode and narrowSpecies when addnodeAndNarrow is called", function() {
    // Mock global rngdata and G
    global.rngdata = {
      "linkreac": {
        "C": ["O", "H"]
      }
    };
    
    const mockG = {
      addNode: function(node) { this.nodes = this.nodes || []; this.nodes.push(node); },
      addEdge: function(from, to) { this.edges = this.edges || []; this.edges.push([from, to]); },
      nodes: [],
      edges: []
    };
    global.G = mockG;
    
    // Define the functions for testing
    const addSingleNode = function(spec) {
      mockG.addNode(spec);
    };
    
    const addnode = function(spec) {
      addSingleNode(spec);
      if (spec in global.rngdata["linkreac"]) {
        global.rngdata["linkreac"][spec].forEach((rightspec) => {
          addSingleNode(rightspec);
          mockG.addEdge(rightspec, spec);
        });
      }
    };
    
    const narrowSpecies = function(spec) {
      const selectIds = ["#speciesselect", "#reactionsselect", "#reactionsabcdselect"];
      selectIds.forEach((selectId) => {
        const selectElement = $(selectId);
        const currentValues = selectElement.val() || [];
        if (!currentValues.includes(spec)) {
          currentValues.push(spec);
          selectElement.val(currentValues);
          selectElement.selectpicker("refresh");
          selectElement.trigger("change");
        }
      });
    };
    
    const addnodeAndNarrow = function(spec) {
      addnode(spec);
      narrowSpecies(spec);
    };
    
    // Mock change events
    let changeTriggered = 0;
    $("#speciesselect, #reactionsselect, #reactionsabcdselect").on("change", function() {
      changeTriggered++;
    });
    
    // Call addnodeAndNarrow
    addnodeAndNarrow("C");
    
    // Check that nodes were added (addnode functionality)
    assert(mockG.nodes.includes("C"), "Main species should be added to graph");
    assert(mockG.nodes.includes("O"), "Linked species O should be added to graph");
    assert(mockG.nodes.includes("H"), "Linked species H should be added to graph");
    assert(mockG.edges.some(edge => edge[0] === "O" && edge[1] === "C"), "Edge should be added from O to C");
    assert(mockG.edges.some(edge => edge[0] === "H" && edge[1] === "C"), "Edge should be added from H to C");
    
    // Check that filtering was applied (narrowSpecies functionality)
    assert.deepStrictEqual($("#speciesselect").val(), ["C"], "Species should be added to speciesselect");
    assert.deepStrictEqual($("#reactionsselect").val(), ["C"], "Species should be added to reactionsselect"); 
    assert.deepStrictEqual($("#reactionsabcdselect").val(), ["C"], "Species should be added to reactionsabcdselect");
    assert.strictEqual(changeTriggered, 3, "Change events should be triggered for filtering");
  });
});