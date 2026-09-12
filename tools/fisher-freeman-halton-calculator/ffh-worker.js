'use strict';
let importError = null;
try { importScripts('/assets/js/ffh-core.js'); } catch (error) { importError = error.message; }
self.onmessage = function (event) {
  const { runId, table, alphaText } = event.data;
  try {
    if (importError) throw new Error('The calculation engine could not be loaded. ' + importError);
    if (!self.FFH || typeof self.FFH.calculate !== 'function') throw new Error('The calculation engine is unavailable.');
    const result = self.FFH.calculate(table, { alphaText });
    self.postMessage({ runId, result });
  } catch (error) {
    self.postMessage({ runId, result: { complete: false, status: 'error', error: error.message } });
  }
};
