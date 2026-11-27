# GUI glossary

Terms used in the pancake GUI for bearing engineers. One row per term: the UI
term, its definition, symbol and unit, and the term used in the final
presentation deck (Korean), for the user guide.

Angles: the fixed frame has x to the right and y up. Fixed-frame angles grow
counterclockwise, which is the direction of rotation for a positive speed.

| UI term | Definition | Symbol, unit | Deck term |
|---|---|---|---|
| Bearing diameter | Nominal bore (shaft) diameter | D = 2R, mm | 내경 |
| Bearing length | Axial length of the film | L, mm | 길이 |
| Radial clearance | Bore radius minus shaft radius | c, µm | 반경 틈새 |
| Oil feed (hole / groove) | Where oil enters the film: a hole of given diameter, or an axial groove over the full length | – | 급유 구멍 / 급유 그루브 |
| Feed angle | Centre of the feed, fixed frame | deg | 급유부 위치 |
| Feed supply overpressure | Feed pressure above the ambient pressure | bar | 급유 압력 |
| Tilt angle, tilt direction | Angle of the shaft axis to the bearing axis, and the direction it leans to | γ (mrad), θ_γ (deg) | 기울어짐 γ |
| Speed | Shaft rotational speed | N, rpm | 회전 속도 |
| Load, load direction | Force on the shaft and the fixed-frame direction it acts in (270° = downward) | W, N | 하중 |
| Specific load | Load per projected area | P = W/(LD), bar | 면압 |
| Ambient pressure | Pressure around the bearing (chamber, shell); both bearing ends are open to it | p_amb, bar | 주변 압력 (챔버 압력) |
| Ambient temperature | Temperature of the surroundings and of the oil entering the film | T_amb, °C | 주변 온도 |
| Moving mass | Mass that moves with the bearing film in a transient case | kg | – |
| Load table | Load (magnitude and direction, or x and y) against crank angle or time | – | 하중 데이터 |
| Oil-refrigerant pair | Oil and dissolved refrigerant with their property table | – | 오일-냉매 |
| Refrigerant content | Mass fraction of refrigerant dissolved in the oil | w, wt% | 냉매 함량 |
| Saturated at ambient pressure | Content at which the oil is in equilibrium with refrigerant at p_amb and T_amb | – | 포화 오일 |
| Saturation pressure | Pressure below which dissolved refrigerant comes out of solution, at the given content and temperature | p_sat, bar | 포화압력 |
| Pressure margin | Ambient minus saturation pressure | p_amb − p_sat, bar | 압력 여유 |
| Pressure drop ratio | How far the gas-free film pressure falls below ambient, per specific load | Φ = (p_amb − p_min)/P | 압력 강하 Φ |
| Onset chart | Φ against S from gas-free runs of the tester bearing; gas release starts when (p_amb − p_sat)/P falls below Φ(S) | – | GC 시작 조건 지도 |
| Gas release expected / borderline / not expected | Onset verdict: margin ratio below the curve / within ±20 % of it / above it. Expected is a design risk (red) | – | GC 시작 / 기체 방출 없음 |
| Sommerfeld number | Bearing duty parameter | S = μ N L D (R/c)² / W, N in rev/s | Sommerfeld 수 S |
| Vapour cavitation (VC) | The oil itself boils at its vapour pressure | VC | VC (증기 캐비테이션) |
| Gaseous cavitation (GC) | Dissolved refrigerant comes out of solution where the film pressure falls below p_sat | GC | GC (기체 캐비테이션) |
| Dissolved-gas release | Model setting. On: GC+VC model, the cavitated film is held at p_sat. Off: VC-only model, the refrigerant stays dissolved | – | GC+VC 모델 / VC 모델 |
| Gas model: pressure plateau / gas transport | Advanced. Plateau tracks the cavitated region: the cavitated film is held at p_sat. Transport also tracks how much refrigerant is released: oil, dissolved and free refrigerant are conserved and solved together. Different goals, not a shortcut and a full version. Transport is used automatically for thermal cases with refrigerant | – | – |
| Cavitated area | Fraction of the film area whose film content is below 1 (filled partly by gas or vapour). The deck metric: 57 % for the compressor example | % | 버블 영역 (유막 면적 대비) |
| Gas fraction | Volume fraction of gas or vapour in the gap; its 0.5 contour outlines the cavitated region in the maps | α | 기체 부피분율 α |
| Minimum film thickness | Smallest gap | h_min, µm | 최소유막두께 |
| Eccentricity ratio | Shaft-centre offset per radial clearance | ε = e/c | 편심률 |
| Attitude angle | Angle from the load line to the line of centres (journal offset relative to the bearing) | φ, deg | attitude angle φ |
| Angle from maximum film | Angle measured from the maximum film thickness in the direction of rotation | θ, deg | θ (최대 유막두께 기준) |
| Angle from feed | Angle measured from the feed centre in the direction of rotation | θ_in, deg | θ_in (급유부 기준) |
| Angle zero | Chosen zero of every angle axis: fixed frame, maximum film, minimum film, feed, opposite the feed | – | – |
| Friction torque | Viscous torque on the shaft | M, N m | 유막 토크 / 마찰 토크 |
| Friction power loss | Torque times angular speed | P_f = Mω, W | 마찰 손실 |
| Peak pressure | Highest film pressure | p_max, bar | 최고 압력 |
| Shear share of film force | Viscous part of the total film force | F_v/F, % | – |
| Surface roughness | RMS roughness of shaft and bore; used only for Λ (no contact model) | Rq, µm | 표면 거칠기 σ |
| Film thickness ratio | Minimum film per combined roughness; below 3 mixed lubrication is likely (red) | Λ = h_min/σ, σ = (Rq₁² + Rq₂²)^½ | 유막두께비 Λ |
| Load balance | Iteration that moves the shaft until the film force carries the load | – | 하중 평형 |
| Iterative solution | Repeated solving until the changes fall below the tolerances (Advanced: FBNS, Newton, Picard, Krylov) | – | 반복 해석 |
| Tolerance | Size of the remaining change or imbalance at which an iteration stops | – | 수렴 판정 |
| Mesh level | Cells around and along the film: coarse 180 × 48, medium 360 × 96, fine 720 × 192 (the deck grid) | n_θ × n_z | 격자 |
| Thermal analysis | Film, housing and shaft temperatures solved (steady only) | – | 열 해석 |
| Steady / transient | Fixed operating point / load varying in time | – | 정상 / 과도 해석 |
| Comparison mode | Runs the case with and without gas release and shows the differences in % | – | GC+VC 대비 VC |
| Boundary zones | Bearing ends, oil feed, shaft surface, bearing surface, housing outer surface | – | 경계 조건 |
| Status light | Green: complete. Yellow: input needed. Red: conflict; the case will not run | – | – |
| Copy command | The exact terminal command that reproduces the case | – | – |
