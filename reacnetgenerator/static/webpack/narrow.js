// SPDX-License-Identifier: LGPL-3.0-or-later

const {getFormula} = require("./formula.js");

const selectIds = [
  "#speciesselect",
  "#reactionsselect",
  "#reactionsabcdselect",
];

/**
 * Add an option for a species that was omitted from the capped select list.
 *
 * Report generation intentionally renders at most 65,536 options, while
 * reaction results can still contain species beyond that limit. A clicked
 * species must exist as an option before jQuery's `val` can select it.
 *
 * @param {HTMLSelectElement} select native select element
 * @param {string} spec species identifier/SMILES
 */
function ensureSpeciesOption(select, spec) {
  const hasOption = Array.from(select.options)
                        .some(
                            (option) => option.value === spec,
                        );
  if (hasOption)
    return;

  const option = select.ownerDocument.createElement("option");
  option.value = spec;
  const formula = getFormula(spec);
  option.textContent = formula ? `${formula} ${spec}` : spec;
  select.appendChild(option);
}

/**
 * Narrow all report sections to include a clicked species.
 *
 * Existing selections are preserved. The report installs jQuery on the
 * global object before this function can be called from the rendered page.
 *
 * @param {string} spec species identifier/SMILES
 */
function narrowSpecies(spec) {
  selectIds.forEach((selectId) => {
    const selectElement = globalThis.$(selectId);
    const select = selectElement[0];
    if (!select)
      return;

    const currentValues = selectElement.val() || [];
    if (currentValues.includes(spec))
      return;

    ensureSpeciesOption(select, spec);
    selectElement.val(currentValues.concat(spec));
    selectElement.selectpicker("refresh");
    selectElement.trigger("change");
  });
}

module.exports = {narrowSpecies};
