# 물성·캐비테이션 모델 개편 정리 (작년 24년 → 올해 25년)

> 작성 목적: 주간 보고 / 물성치 섹션 초안. 각 항목을 **작년식 → 올해식 → 선정 근거 → 타당성(솔직) → 미구현** 구조로 정리.
> 데이터 출처 확인: `01. 냉매 정보 / 점도그래프_값추정_PTSV_PZ68S_R290_V2.0_230703.xlsm` (박상백, LG전자), pancake `config_r290_pz68.txt` / `PHYSICS.md` / `AUDIT.md`.

---

## 0. 모델 골격 — 혼합물을 몇 상(phase)으로 보는가 (이번에 확정)

| 그림 | 액체상 | 냉매 취급 | 비고 |
|---|---|---|---|
| 작년(교수님) | 오일 | 용존기체(dissolved) + 방출기체(free) | 용해된 냉매를 **기체**로 취급 |
| 회사 | 오일 + 액냉매 | 용해=액상, 과포화 시 기화→버블 | 용해된 냉매를 **액상**으로 취급 |
| **확정** | **오일 + 액냉매(용액)** | **액냉매(용해) + 기냉매(버블) 3성분** | 용존 냉매를 액상으로 일관 |

**선정 근거.** PTSV 선도는 "오일에 냉매가 녹은 **용액**의 점도·용해도를 직접 측정"한 데이터다. 즉 용해된 냉매의 부피 기여가 이미 액상으로 반영돼 있다. 열역학적으로도 *액체에 녹은 기체 분자의 부분몰부피는 기체상이 아니라 액상에 가깝다.* 따라서 "용해 = 액상"이 측정과 일치하며, 작년의 "용존냉매에 이상기체 밀도" 처리가 실제 오류였다.

**결과 정의.**
- 액체 용액 = 오일 + 용해 냉매, 물성 $\rho,\mu$ 는 (가능하면) PTSV 차트 직접 사용.
- 기체상 = 국소 농도가 포화용해도 $w_{sat}(p,T)$ 를 넘을 때 방출된 냉매 버블. 별도 추적, 두께방향 부피분율 $\alpha$ 로 유효물성에 반영.
- 오일 자체 증기(vaporous)는 부차적 → 기본 끔.

---

## 1. Cavitation 모델: Kunz → Elrod–Adams (JFO)

**작년.** Kunz vaporous 모델. 부피분율을 직접 수송·조작 → **질량 비보존**. 검증 수단 없음.

**올해.** Elrod–Adams JFO. 단일 변수 film fraction $\theta$ 로 full-film과 cavitated를 통합.

$$\rho = \theta\,\rho_{liq},\qquad
\begin{cases}\theta \ge 1,\ p > p_{cav} & \text{(full film)}\\[2pt]\theta < 1,\ p = p_{cav} & \text{(cavitated)}\end{cases}$$

스위치 함수 $g$ 로 Poiseuille 항을 cavitated 영역에서 꺼서, 한 식으로 rupture–reformation을 자동 처리:

$$\frac{\partial}{\partial x}\!\Big(\frac{g\,\rho h^3}{12\mu}\frac{\partial p}{\partial x}\Big)+\frac{\partial}{\partial z}\!\Big(\frac{g\,\rho h^3}{12\mu}\frac{\partial p}{\partial z}\Big)=\frac{U}{2}\frac{\partial(\rho h)}{\partial x}+\frac{\partial(\rho h)}{\partial t}$$

**근거.** Elrod–Adams는 저널베어링 cavitation의 **질량보존 표준**(JFO 이론의 이산형). 작년 Kunz가 위배하던 유막 질량보존을 만족.

**타당성(확인됨).** 매 스텝 액체/기체 전역 질량수지 잔차 $\sim 10^{-12}$, 격자수렴 확인(cav 면적 5.17→4.83→4.80 %, 60×20/120×40/240×80), 테스트 19/19. → "버블 위치를 그리던 코드 → 질량보존이 입증되는 예측 코드".

---

## 2. 포화용해도 $w_{sat}(p,T)$

**작년.** FW68D / **R410A** PTSV 선도 추출점을 압력·온도 멱함수로 fitting.

$$w_{sat}(p,T)=a\,p^{\,b},\quad a=1060\,T^{-1.389},\ \ b=0.927-5.652\times10^{-9}(T-35)^9$$

**올해.** **R290 / PZ68S** 차트의 (40 °C, 1.0 MPa) 기준점에 Henry + van't Hoff.

$$c_{sat}(p,T)=H(T)\,p,\qquad H(T)=H_{ref}\exp\!\Big[E_H\Big(\tfrac1T-\tfrac1{T_{ref}}\Big)\Big]$$

(코드: `H_ref = 2.175\times10^{-7}\,\text{Pa}^{-1}`, `E_H = 1800\,K`. → 1.0 MPa에서 $c_{sat}=0.2175$ = 21.75 %.)

**근거.** Henry 법칙은 **희박 기체–액체 용해 평형의 표준**(농도 ∝ 분압), van't Hoff 항이 Henry 상수에 열역학적 온도의존성 부여. 무엇보다 **대상 냉매 R410A→R290 교체**가 핵심.

**타당성(절반).** 냉매 교체는 명확한 개선. 그러나 선형 Henry는 **1점 접선**이라 고압부 미검증 — 4 MPa로 외삽하면 $H\!\cdot\!p \approx 0.87$ (87 %, 비물리). 실제 차트는 곡률(저압 과소·고압 포화)을 가지므로 이 구간이 결정적. → **TABLE 필요(미구현, §8)**.

---

## 3. 점도 $\mu$

**작년.** 40 / 100 °C 두 점을 지나는 지수함수. 40 °C 이하는 과대평가 방지를 위해 **기울기를 연속($C^1$)으로 맞춘 선형식**으로 접합(꺾임 아님). 혼합은 **부피분율** 가중.

$$\mu_{oil}(T)=A\,e^{-B T}\ (T>40^\circ\text{C}),\qquad \mu=\text{linear}(T)\ \text{with matched slope}\ (T\le40^\circ\text{C})$$

한계: 압력 의존성 없음, 용존냉매 희석효과 없음.

**올해(액상).** 점도를 (온도·용존농도·압력) 세 표준효과의 곱으로.

$$\mu_l(T,c_d,p)=\mu_{ref}\,\underbrace{\exp\!\Big[E_\mu\big(\tfrac1T-\tfrac1{T_{ref}}\big)\Big]}_{\text{Andrade (온도)}}\ \underbrace{\exp(a_c\,c_d)}_{\text{용존냉매 희석}}\ \underbrace{\exp[\alpha_p(p-p_{ref})]}_{\text{Barus (압력)}}$$

(코드: $E_\mu=4000\,K$, $a_c=-11.40$, $\alpha_p=1.5\times10^{-8}\,\text{Pa}^{-1}$. → 40 °C, 1.0 MPa에서 $\nu=6.843$ cSt.)

**올해(2상 혼합).** 자유기체가 생기면 유효점도를 **질량분율(quality) 기반 McAdams**로(작년 부피분율 → 질량분율 전환).

$$\frac{1}{\mu_{mix}}=\frac{x}{\mu_g}+\frac{1-x}{\mu_l},\qquad x=\frac{m_g}{m_g+\rho_l\theta h}\ \text{(질량분율)}$$

(Krieger–Dougherty 패킹 보정 포함, $\mu_g\approx8.7\times10^{-6}$ Pa·s = 포화 프로판 증기.)

**근거(항목별).**
- **Andrade $\exp(E_\mu/T)$** — 액체 점도–온도의 표준식. 작년의 지수+선형 패치를 단일 물리식으로 대체(저속·저온 운전 목표에서 신뢰성).
- **희석 $\exp(a_c c_d)$** — 작년 점도식은 용존농도와 무관. 측정된 점도강하(순수오일 대비)를 재현하기 위해 도입.
- **McAdams 질량분율** — 균질 2상 점도의 표준. 작년 부피분율 가중은 **고 void에서 발산/비물리**, 질량분율은 두 단상 극한($\mu_l,\mu_g$)을 모두 만족 → 이전 모델이 못 잡던 고 void 영역을 위해 선정.

**타당성.** Andrade·희석·질량분율 전환은 근거 명확·정량 확인. **단 Barus $\alpha_p$ 는 PZ68 측정값이 아닌 광유 일반값** → "근거 있는 디폴트"일 뿐 미검증(§8).

---

## 4. 밀도 $\rho$

**작년.** 오일 880 kg/m³ 상수, **용존냉매 = 오일밀도**, 방출기체 = 이상기체, 혼합 **부피분율** 가중.

**올해(용액).** 용존냉매를 액상 분자부피로 보고 질량–부피 혼합.

$$\frac{1}{\rho_{sol}}=\frac{1-c_g}{\rho_{oil}}+\frac{c_g}{\rho_{g,liq}},\qquad \rho_{g,liq}\approx500\ \text{kg/m}^3$$

**올해(2상 혼합, cavitation).** 버블이 생기면 혼합밀도가 지배적으로 변하며, 이는 Elrod $\theta$ 가 담는다.

$$\rho_{mix}=(1-\alpha)\rho_l+\alpha\,\rho_g\;\;\Longleftrightarrow\;\;\rho=\theta\rho_{liq}$$

**근거.** 용존냉매는 액상에 분자 단위로 녹아 있으므로 액상 비체적으로 가산(Amagat형 부피가산 = 액–액 혼합 표준). 작년의 이상기체 밀도 오류 제거.

**타당성(솔직).** *용액* 밀도 보정은 개념적으로 옳으나 Reynolds 해에 주는 영향은 점도·cavitation보다 작음 → "틀린 가정의 정합성 보정". 반면 *2상(버블) 밀도*는 무시 불가 — **cavitation 모델 본체**(프로판 증기/액체 밀도비 ~0.01–0.04라 $\alpha$ 조금에도 밀도 1–2자릿수 하락). PTSV 차트에는 밀도가 없으므로 밀도는 처음부터 끝까지 모델 담당.

---

## 5. 혼합 비율: 부피분율 → 질량분율

작년은 점도·밀도 혼합을 **부피분율(volume fraction)**로 계산. 올해는 **질량분율(quality)**로 전환. 질량분율 혼합이 두 단상 극한을 만족하고 고 void에서 물리적이라는 점이 §3·§4의 공통 근거.

---

## 6. 기체 방출이 점도·밀도에 주는 효과 — 무시 가능한가?

차트(오일+액냉매)는 **액체 한 상**만 준다. 기체 방출은 그 위에 얹히는 **별개의 2상 효과**라 차트로는 안 나온다.

- 기체 방출 시 점도가 **두 방향**으로 바뀜:
  - (a) 액상 점도 ↑ : 희석제(용존냉매)가 빠져 $c_d$ ↓ → $\nu(T,c_d)$ 가 순수오일 쪽으로 회복 (차트가 $c_d$ 추적으로 자동 처리).
  - (b) 혼합 점도 ↓ : 기포가 섞여 유효점도가 기체 쪽으로 하락 (차트 밖, McAdams가 처리).
- **밀도**: 무시 불가 = cavitation 본체(§4).
- **점도(2상)**: Reynolds의 Poiseuille 항 $\rho h^3/12\mu$ 는 cavitated 영역에서 $p\approx$ 일정이라 거의 0 → **압력/하중엔 2차**. 반면 마찰응력 $\tau=\mu U/h$ 엔 그대로 들어가 → **마찰토크·rupture 경계엔 1차**. ($\mu_g/\mu_l\sim1.6\times10^{-3}$.)

결론: 액상을 차트로 정확히 알수록 2상 혼합의 "액체 끝점"이 정확해져 기체효과 예측도 같이 좋아진다. 둘 다 필요하되 비중이 다르다.

---

## 7. 물성 외 — 작년 골격 대비 (참고)

| 항목 | 작년 | 올해 | 근거 |
|---|---|---|---|
| 선형 솔버 | Jacobi | PETSc Krylov | 수렴·견고성 |
| 대류 이산화 | upwind 고정 | upwind / TVD / type-diff 선택 | front 가짜확산(수렴차수 ~1.09→~1.99) |
| 검증 | 격자수렴(1 %) | 격자수렴 + 매 스텝 질량수지(잔차 1e-12) | 신뢰성 입증 |
| Transient | (capacity term 누락) | $\rho_0 h\,\partial\theta/\partial t$ 저장항 복원 | 과도해석 정확성 |
| 인터페이스 | Linux config.txt | Windows GUI (ImGui/ImPlot) | 사용자 매뉴얼 대상 |

지배방정식(Reynolds 윤활 + 에너지 보존), 저널 다이나믹스, MPI 병렬은 작년 골격 유지·정비.

---

## 8. 아직 구현 안 된 것 (미구현 목록)

> `config_r290_pz68.txt`에서 직접 확인. 세 물성 모두 TABLE 옵션은 존재하나 표가 비어 있음:
> `oil_gas_solution_model = HENRY` (옵션에 TABLE), `viscosity_model = EMPIRICAL_CORRELATION` (옵션에 TABLE), `density_model = MASS_VOLUME_MIXING` (옵션에 TABLE),
> `solubility_table = ` (빈칸), `density_table = ` (빈칸), `viscosity_table = ` (빈칸).
> 주석: *"Henry + exp(a_c c_d) viscosity are local tangents calibrated at (10 bar, 40 C) only; measured R290/PZ68 isotherm table data is pending."*

| # | 미구현 항목 | 현재 상태 | 해야 할 일 | 우선도 |
|---|---|---|---|---|
| 1 | **PTSV TABLE 입력** | 표 3종 비어 있음 (1점 접선으로 동작) | xlsm에서 $w_{sat}(p)$, $\nu(c_d)$ 추출 → `*_table` 채움 + 현재 모델 대비 편차맵 | **최우선** |
| 2 | **TABLE 차원** | 코드 TABLE은 1D (용해도=압력키, 점도/밀도=용존질량분율키) | 차트는 2D $(p,T)$ → 기준 등온선 채택 or 2D 확장 결정 | 높음 |
| 3 | **방출기체 압력 커플링** (AUDIT WP-1/2) | released gas가 압력 해석에 **부피적으로 비활성**, void onset이 고정 $p_{cav}$ | $p_{sat}(T,c_d)$ 기반 onset + 연속방정식 커플링 (문헌상 하중 −3~+21 %) | **최우선** |
| 4 | **용존=액상 일관화 마무리** | `dissolved_gas_liquid_density=500`로 절반만 됨 | 용액 $\rho,\mu$ 출처를 mixing-fit → PTSV로 전환(§0 확정 반영) | 높음 |
| 5 | **Barus $\alpha_p$ 검증** | `1.5e-8` = 광유 일반값 | PZ68 측정 or 문헌 보정, 아니면 영향 평가 후 명시 | 중 |
| 6 | **$\dot m$ 방출률 상수** | `gas_mass_transfer_rate = 500` (임의) | 민감도 스터디 후 문헌/실험 보정 | 중 |
| 7 | **오일 vaporous 처리** | 기본 끔(확정) | 옵션 유지 여부만 결정 | 낮음 |
| 8 | **시험기 가시화 대조** | 없음 | cavitation 영역 실험 비교(검증 최종 단계) | 중장기 |

---

## 9. 한 줄 종합 (솔직)

- **확실한 개선(물리 근거 명확):** 냉매 R410A→R290, Kunz→Elrod JFO(질량보존), 용존냉매 점도 희석, 부피분율→질량분율(고 void 발산 제거), 이상기체 밀도 오류 제거.
- **방향만 맞고 미검증:** 용해도 1점 선형 Henry(고압부), Barus 압력항, 방출률 상수.
- **완결 조건:** PTSV 차트 전체를 TABLE로 입력하고(§8-1,2), 방출기체를 압력에 커플링(§8-3)해야 "물성·캐비테이션 모델 개편"이 완결된다.

### 데이터 커버리지 (PTSV xlsm 확인 결과)
- 신뢰범위: 압력 **0.1–4 MPa**, 온도 **−60 ~ +160 °C**, 점도 ~1000 mm²/s, R²=0.9995.
- 출력: (P,T) → 점도 + 냉매용해도. **밀도는 없음.**
- 성질: 차트 이미지를 픽셀 디지타이즈 → 다항식 피팅한 **추정 모델**(11개 등압선, 각 ~10–20점). 측정 실측이 아니라 *차트 1장의 재구성* — 1점 접선보다는 압도적이나 원본 차트 정확도를 물려받음.
- 포화면(saturation surface)을 주지만, 등온선을 따라 $(\nu_{sat},w_{sat})$ 를 교차플롯하면 불포화 $\nu(T,c_d)$ 까지 복원 가능 → 코드의 `viscosity_table`(용존질량분율 키)과 직접 호환.
