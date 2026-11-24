import { reverseComplement } from "../dna/alphabet.js";
import { Rng } from "./random.js";

export interface GenomeOptions {
  length: number;
  gc?: number;
  repeatFamilies?: number;
  repeatLength?: number;
  repeatCopies?: number;
  repeatDivergence?: number;
  seed?: number;
}

export interface PlantedRepeat {
  family: number;
  start: number;
  length: number;
  reverse: boolean;
}

export interface SimulatedGenome {
  sequence: string;
  repeats: PlantedRepeat[];
}

function randomBase(rng: Rng, gc: number): number {
  const r = rng.next();
  if (r < gc / 2) return 67;
  if (r < gc) return 71;
  if (r < gc + (1 - gc) / 2) return 65;
  return 84;
}

const OTHER_BASES: Record<number, number[]> = {
  65: [67, 71, 84],
  67: [65, 71, 84],
  71: [65, 67, 84],
  84: [65, 67, 71]
};

export function mutateBase(rng: Rng, base: number): number {
  const options = OTHER_BASES[base] ?? OTHER_BASES[65];
  return options[rng.int(3)];
}

export function simulateGenome(options: GenomeOptions): SimulatedGenome {
  const length = Math.floor(options.length);
  const gc = options.gc ?? 0.5;
  const families = options.repeatFamilies ?? 0;
  const repeatLength = options.repeatLength ?? 0;
  const copies = options.repeatCopies ?? 2;
  const divergence = options.repeatDivergence ?? 0;
  if (length < 1) throw new Error("genome length must be positive");
  if (!(gc >= 0 && gc <= 1)) throw new Error("gc must be between 0 and 1");
  if (families > 0 && (repeatLength < 1 || copies < 2)) throw new Error("repeats need a positive length and at least 2 copies");
  if (families * copies * repeatLength > length / 2) throw new Error("repeats would cover more than half the genome");
  const rng = new Rng(options.seed ?? 1);
  const bases = Buffer.allocUnsafe(length);
  for (let i = 0; i < length; i++) bases[i] = randomBase(rng, gc);
  const repeats: PlantedRepeat[] = [];
  const occupied: Array<[number, number]> = [];
  const overlaps = (start: number, end: number): boolean => occupied.some(([s, e]) => start < e + 1 && s < end + 1);
  for (let family = 0; family < families; family++) {
    const unit = Buffer.allocUnsafe(repeatLength);
    for (let i = 0; i < repeatLength; i++) unit[i] = randomBase(rng, gc);
    const unitForward = unit.toString("latin1");
    const unitReverse = reverseComplement(unitForward);
    for (let copy = 0; copy < copies; copy++) {
      let placed = false;
      for (let attempt = 0; attempt < 1000 && !placed; attempt++) {
        const start = rng.int(length - repeatLength + 1);
        if (overlaps(start, start + repeatLength)) continue;
        const reverse = copy > 0 && rng.next() < 0.5;
        const source = reverse ? unitReverse : unitForward;
        for (let i = 0; i < repeatLength; i++) {
          let base = source.charCodeAt(i);
          if (divergence > 0 && rng.next() < divergence) base = mutateBase(rng, base);
          bases[start + i] = base;
        }
        occupied.push([start, start + repeatLength]);
        repeats.push({ family, start, length: repeatLength, reverse });
        placed = true;
      }
      if (!placed) throw new Error("could not place repeat copies without overlap; use fewer or shorter repeats");
    }
  }
  repeats.sort((x, y) => x.start - y.start);
  return { sequence: bases.toString("latin1"), repeats };
}
