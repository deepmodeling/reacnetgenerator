/** This module is proposed to select elements in parallel. */

/**
 * Select from a list.
 * Reference: https://segmentfault.com/a/1190000022360080
 * @param {array} select The list to select.
 * @param {array} from The list to select from.
 * @param {function} func The function to check if an element should be selected, which
 *      has two parameters: element, select.
 * @return {array} The selected list.
 */
const search = async (select, from, func) => {
    // if there is no select, return all of list
    if (!select.length) return from;
    const results = await Promise.all(from.map(element => func(element, select)));
    return from.filter((_v, index) => results[index]);
}

const searchspecies = (select, from) => await search(select, from, async(element, select) => {
    return select.includes(element);
});

const searchreaction = (select, from) => await search(select, from, async(element, select) => {
    return [element["l"], element["r"]].some(spec => select.includes(spec));
});

const searchreactionabcd = (select, from) => await search(select, from, async(element, select) => {
    return element["l"].concat(element["r"]).some(spec => select.includes(spec));
});

module.exports = { searchspecies, searchreaction, searchreactionabcd };