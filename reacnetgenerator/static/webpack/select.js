// SPDX-License-Identifier: LGPL-3.0-or-later
/**
 * Select from a list.
 * @param {array} select The list to select.
 * @param {array} from The list to select from.
 * @param {function} func The function to check if an element should be
 *     selected, which
 *      has two parameters: element, select.
 * @return {array} The selected list.
 */
const search =
    (select, from, func) => {
      // if there is no select, return all of list
      if (!select.length)
        return from;
      return from.filter(element => func(element, select));
    }

const searchspecies = (select, from) =>
    search(select, from,
           (element, select) => { return select.includes(element["s"]); });

const searchreaction = (select,
                        from) => search(select, from, (element, select) => {
  return element["l"].concat(element["r"]).some(spec => select.includes(spec));
});

module.exports = {
  searchspecies,
  searchreaction
};
