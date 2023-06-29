// SPDX-License-Identifier: LGPL-3.0-or-later
const assert = require('assert');
const {getFormula} = require('../formula');

/** test getFormula */
describe('Formula', function() {
  describe('getFormula', function() {
    it('uppercase',
       function() { assert.equal(getFormula('[C][C]([H])[C][C][C]'), 'C5H'); });
    it('lowercase', function() {
      assert.equal(getFormula('[H][c][c][c][c][c][c]'), 'C6H');
    });
    it('capital', function() {
      assert.equal(getFormula('[Ca][c][C][c][Ca][c]'), 'C4Ca2');
    });
  });
});
