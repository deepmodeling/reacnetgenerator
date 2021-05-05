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
const search = (select, from, func) => {
    // if there is no select, return all of list
    if (!select.length) return from;
    const searchSync = async (select, from, func) => Promise.all(from.map(element => func(element, select)))
        .then((results) => from.filter((_v, index) => results[index]))
    return await searchSync(select, from, func);
}

const searchspecies = (select, from) => search(select, from, async (element, select) => {
    return select.includes(element);
});

const searchreaction = (select, from) => search(select, from, async (element, select) => {
    return [element["l"], element["r"]].some(spec => select.includes(spec));
});

const searchreactionabcd = (select, from) => search(select, from, async (element, select) => {
    return element["l"].concat(element["r"]).some(spec => select.includes(spec));
});

module.exports = { searchspecies, searchreaction, searchreactionabcd };