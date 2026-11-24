export const BASES = "ACGT";

export const BASE_CODE: Int8Array = (() => {
  const table = new Int8Array(256).fill(-1);
  const pairs: Array<[string, number]> = [
    ["A", 0],
    ["C", 1],
    ["G", 2],
    ["T", 3]
  ];
  for (const [base, code] of pairs) {
    table[base.charCodeAt(0)] = code;
    table[base.toLowerCase().charCodeAt(0)] = code;
  }
  return table;
})();

const COMPLEMENT: Uint8Array = (() => {
  const table = new Uint8Array(256);
  for (let i = 0; i < 256; i++) table[i] = 78;
  const pairs = ["AT", "CG", "GC", "TA", "at", "cg", "gc", "ta", "NN", "nn"];
  for (const pair of pairs) table[pair.charCodeAt(0)] = pair.charCodeAt(1);
  return table;
})();

export function reverseComplement(seq: string): string {
  const n = seq.length;
  const out = Buffer.allocUnsafe(n);
  for (let i = 0; i < n; i++) out[n - 1 - i] = COMPLEMENT[seq.charCodeAt(i) & 0xff];
  return out.toString("latin1");
}

export function reverseString(value: string): string {
  const n = value.length;
  const out = Buffer.allocUnsafe(n);
  for (let i = 0; i < n; i++) out[n - 1 - i] = value.charCodeAt(i) & 0xff;
  return out.toString("latin1");
}

export interface BaseComposition {
  a: number;
  c: number;
  g: number;
  t: number;
  n: number;
}

export function emptyComposition(): BaseComposition {
  return { a: 0, c: 0, g: 0, t: 0, n: 0 };
}

export function addComposition(seq: string, into: BaseComposition): BaseComposition {
  for (let i = 0; i < seq.length; i++) {
    switch (BASE_CODE[seq.charCodeAt(i) & 0xff]) {
      case 0:
        into.a++;
        break;
      case 1:
        into.c++;
        break;
      case 2:
        into.g++;
        break;
      case 3:
        into.t++;
        break;
      default:
        into.n++;
    }
  }
  return into;
}

export function gcFraction(composition: BaseComposition): number {
  const acgt = composition.a + composition.c + composition.g + composition.t;
  return acgt === 0 ? 0 : (composition.g + composition.c) / acgt;
}
