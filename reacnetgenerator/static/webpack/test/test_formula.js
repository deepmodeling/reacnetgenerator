const assert = require('assert');
const { getFormula } = require('../formula');

/** test getFormula */
describe('Formula', function () {
  describe('getFormula', function () {
    it('uppercase', function () {
      assert.equal(getFormula('[C][C]([H])[C][C][C]'), 'C<sub>5</sub>H');
    });
    it('lowercase', function () {
        assert.equal(getFormula('[H][c][c][c][c][c][c]'), 'C<sub>6</sub>H');
    });
    it('capital', function () {
        assert.equal(getFormula('[Ca][c][C][c][Ca][c]'), 'C<sub>4</sub>Ca<sub>2</sub>');
    });
  });
});
