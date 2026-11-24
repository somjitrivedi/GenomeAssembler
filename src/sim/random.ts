export class Rng {
  private a: number;
  private b: number;
  private c: number;
  private d: number;
  private spare: number | null = null;

  constructor(seed: number) {
    let s = seed >>> 0;
    const next = (): number => {
      s = (s + 0x9e3779b9) >>> 0;
      let z = s;
      z = Math.imul(z ^ (z >>> 16), 0x85ebca6b);
      z = Math.imul(z ^ (z >>> 13), 0xc2b2ae35);
      return (z ^ (z >>> 16)) >>> 0;
    };
    this.a = next();
    this.b = next();
    this.c = next();
    this.d = next() | 1;
  }

  nextUint32(): number {
    const result = Math.imul(((Math.imul(this.b, 5) << 7) | (Math.imul(this.b, 5) >>> 25)), 9) >>> 0;
    const t = this.b << 9;
    this.c ^= this.a;
    this.d ^= this.b;
    this.b ^= this.c;
    this.a ^= this.d;
    this.c ^= t;
    this.d = (this.d << 11) | (this.d >>> 21);
    return result;
  }

  next(): number {
    return this.nextUint32() / 4294967296;
  }

  int(maxExclusive: number): number {
    return Math.floor(this.next() * maxExclusive);
  }

  normal(mean: number, sd: number): number {
    if (this.spare !== null) {
      const value = this.spare;
      this.spare = null;
      return mean + sd * value;
    }
    let u = 0;
    while (u === 0) u = this.next();
    const v = this.next();
    const r = Math.sqrt(-2 * Math.log(u));
    this.spare = r * Math.sin(2 * Math.PI * v);
    return mean + sd * r * Math.cos(2 * Math.PI * v);
  }
}
