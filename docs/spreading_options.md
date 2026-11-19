# Spreading the basalt: do we need to buy a spreader?

Draft, 19 November 2026. Tonnages come from the model (placeholder moisture 8%). Prices are public benchmarks and one listing, not quotes.

## The job
- 12 treated plots of 0.64 ha (7.68 ha), each at one fixed rate: 20, 40 or 80 t/ha dry.
- 390 t wet to spread (358 t dry), plus contingency left in the stockpile.
- One spreading event, at the start of the trial. Nothing to spread again within the 30 months.

## What the AS65 is (Agri-Spread AS65-2K)
- Trailed single-axle belt-floor spreader with spinning discs, 6.6 m3 heaped, 9.2 t carrying capacity, four-point load cells.
- Optional ISOBUS auto rate control (variable rate, GPS, section control) or a simpler land-drive version with no electronics.
- Spreads lime to 12 m (manufacturer figure). Not published: maximum application rate, belt discharge rate, tractor power needed.
- Price: no list price published; dealers quote on request. One used 2022 unit was advertised in Aberdeen at about USD 38,700. A new ISOBUS machine will cost more; get a dealer quote. Shipping to Sao Miguel comes on top.

## Do we need variable rate? No.
Variable-rate technology changes the rate *within* a field from a map. Our design needs a uniform rate *within* each plot and a different rate *between* plots. What we need is:
1. **Accurate, verified rate per plot**: a calibrated spreader plus the 5 catch trays per plot already in the design. Load cells help (you weigh what went on each plot) but are not essential if every plot gets a weighed number of loads.
2. **Clean edges**: lime throws 12 m, our alleys are 8 m. Without border control, rock lands in alleys and neighbouring plots. Border spreading (a border disc, or driving the outer bout inward) matters more than variable rate.
3. **High rates**: 80 t/ha is about 10 times a normal lime dose. At 12 m width and 4 km/h a single pass would need about 6.4 t per minute out of the spreader, which is not realistic for a spinner. Plan on several passes (for example 4 x 20 t/ha), which also evens out the spread. 40 t/ha means about 2 to 3 passes, 20 t/ha 1 to 2.

## Design consequence to fix now
config.json assumes plots of "five 16 m spreading bouts". At a 12 m lime width, 80 m is 6.7 bouts. Either confirm the contractor's real working width and adjust the plot or bout plan, or accept a narrower final bout. The catch-tray layout should follow the bout pattern.

## Loads and time
| Rate | Wet t per plot | 9.2 t loads per plot | Plots |
|---|---|---|---|
| 20 t/ha | 13.9 | 2 | 4 |
| 40 t/ha | 27.8 | 4 | 4 |
| 80 t/ha | 55.7 | 7 | 4 |
| Total | 390 | 52 | 12 |

At about 30 minutes per load (load at the stockpile, drive, spread), that is roughly 26 machine hours: **3 to 4 working days with one machine**, plus a half day for calibration on a test strip with trays.

## Three options
| Option | Cost (benchmark) | For | Against |
|---|---|---|---|
| **A. Hire a contractor with a lime spreader** | Ireland 2026 guide EUR 10.50 to 12.50/t spread (FCI); UK 2025 about GBP 7.47/t (NAAC). On 390 t: **about EUR 4,100 to 4,900**, before any high-rate or multi-pass surcharge. Azores rates unknown | Cheapest, experienced operator, no asset to own | Need a local contractor with border control and willingness to do multi-pass, high-rate work on small plots |
| B. Hire the machine only | Unknown; ask dealers and the local co-op | Control over the operator | Still need a tractor of suitable size, a trained driver, and transport to the island |
| C. Buy an AS65 | Used example about USD 38,700; new ISOBUS more; plus shipping and a tractor | Available for any later phase | 10 times option A or more for 4 days of work; check whether the call even funds equipment purchase or only depreciation over the project |

**Recommendation: option A.** Buying is not justified for one 4-day job on 7.7 ha. Revisit buying only if the commercial phase after the trial needs a spreader regularly; then the AS65's load cells and ISOBUS start to pay off.

## To ask before the application
1. Luis: who spreads lime on the farm now, with what machine, at what width; is there a contractor on Sao Miguel who would take this on?
2. Contractor: working width for damp basalt fines, border control, maximum rate per pass, price per tonne or per hour for multi-pass work, dust control.
3. Tecnovia: bulk density and moisture of the fines (sets how many tonnes fit in a load) and whether the quarry loader can load at the farm stockpile.
4. Juan: is a multi-pass, cross-direction pattern acceptable for the protocol, and how the catch trays should be placed.

## Sources
- Agri-Spread AS65-2K product page: https://www.agrispread.com/products/model/as65-2k
- Used AS65 listing (Machinio): https://www.machinio.com/cat/agri-spread
- FCI contracting charges guide 2026 (Agriland): https://www.agriland.ie/farming-news/rates-in-fci-contracting-charges-guide-up-4-for-2026/
- FCI contracting charges guide 2025: https://www.farmcontractors.ie/post/fci-contracting-charges-guide-2025
- NAAC rates (The Scottish Farmer): https://www.thescottishfarmer.co.uk/news/26099857.slurry-spreading-among-biggest-contracting-costs-year/
