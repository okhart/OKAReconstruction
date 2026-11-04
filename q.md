## Реализация многократного рассеяния в ACTS v9.0.0

### 1. Расчёт среднеквадратичного угла рассеяния

**Файл:** `Core/src/Material/Interactions.cpp`

- `theta0Highland()` — формула Хайленда для большинства частиц.

$$
\theta_0 =
\frac{13.6\,\text{МэВ}}{\beta p}|q|
\sqrt{\frac{x}{X_0}}
\left[
1+0.038\ln\left(\frac{x}{X_0}\frac{q^2}{\beta^2}\right)
\right].
$$

- `theta0RossiGreisen()` — формула Росси — Грайзена для электронов и позитронов.

$$
\boxed{
\theta_0 =
\frac{17.5\,\mathrm{MeV}}{\beta p}|q|
\sqrt{\frac{x}{X_0}}
\left[
1+0.125\log_{10}\left(
10\frac{x}{X_0}
\right)
\right]
}
$$

- `computeMultipleScatteringTheta0()` — выбор формулы и расчёт угла рассеяния по параметрам частицы и материала.

### 2. Получение характеристик материала

**Файл:** `Core/include/Acts/Propagator/detail/PointwiseMaterialInteraction.hpp`

Функция `evaluateMaterialSlab()`:
- Получает свойства материала пересечённой поверхности.
- Вычисляет поправку на угол падения `pathCorrection`.
- Корректирует пройденной толщину материала.

### 3. Вычисление дисперсий

**Файл:** `Core/src/Propagator/detail/PointwiseMaterialInteraction.cpp`

- `evaluatePointwiseMaterialInteraction()` — запускает расчёт эффектов, связанных со свойствами материалов.
- `covarianceContributions()` — вычисляет добавки к ковариации.

### 4. Обновление ковариационной матрицы

**Файл:** `Core/include/Acts/Propagator/detail/PointwiseMaterialInteraction.hpp`

Функции:
- `updateState()` — обновляет состояние частицы и ковариационную матрицу.
- `updateVariance()` — добавляет или вычитает рассчитанную дисперсию.

При прямом распространении:

$$
P_{\phi\phi}'=P_{\phi\phi}+\sigma_\phi^2,
$$

$$
P_{\theta\theta}'=P_{\theta\theta}+\sigma_\theta^2.
$$

(P - ковариационная матрица).


### 5. Использование в CKF

**Файл:** `Core/include/Acts/TrackFinding/CombinatorialKalmanFilter.hpp`

В ACTS траектория анализируется пошагово, поэтому на каждом шаге проводится вычисление ковариационной матрицы на основе свойств веществ и поверхностей.

- `materialInteractor()` — связывает обработку материала с алгоритмом CKF.
- `filter()` — вызывает обработку материала на поверхностях детектора.
