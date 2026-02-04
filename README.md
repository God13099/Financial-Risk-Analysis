# Financial Risk Analysis: Portfolio Volatility & Tail Dependence

### 📊 Project Overview
This project performs a comprehensive risk analysis on a 50/50 equity portfolio consisting of **Santander Bank Polska (SPL)** and **KRUK S.A. (KRU)**. 

Using advanced econometric methods, the analysis moves beyond standard Gaussian assumptions to model **fat tails**, **volatility clustering**, and **non-linear dependence** between the banking and debt management sectors. The project concludes with a rigorous backtesting framework compliant with **Basel II/III regulatory standards**.

### 🛠️ Key Methodologies
The project is implemented in **R** and covers the following quantitative workflows:

1.  **Univariate Modeling**: 
    * Detected heavy tails and volatility clustering in asset returns.
    * Selected **eGARCH(1,1) with Student-t distribution** as the optimal model based on AIC/BIC criteria.
2.  **Dependence Structure**:
    * Modeled joint distribution using **Copula** theory.
    * **t-Copula** was selected over Gaussian/Archimedean families (highest Log-Likelihood), confirming significant tail dependence between assets.
3.  **Risk Forecasting (VaR & ES)**:
    * Estimated **Value at Risk (VaR)** and **Expected Shortfall (ES)** at 1% and 5% levels.
    * Compared Static models (Historical, Normal) vs. Dynamic models (EWMA, GARCH, DCC-GARCH).
4.  **Backtesting & Validation**:
    * **Statistical Tests**: Kupiec (Frequency), Christoffersen (Independence), Pelletier (Duration), and Frey-McNeil (Magnitude).
    * **Regulatory Tests**: Implemented the **Basel Traffic Light** approach (rolling 250-day window).

### 📉 Key Findings
* **Gaussian Failure**: The Normal distribution consistently underestimates tail risk (ES -5.20% vs Historical -7.24%), proving dangerous for stress testing.
* **The "Iron Shield"**: **DCC-GARCH** was the only model to remain in the **Basel Green Zone** (0-4 exceedances) throughout the 8-year history.
* **Regime Shifts**: Simple models like **EWMA** failed during market crashes (e.g., COVID-19), falling into the "Red Zone" and incurring maximum regulatory capital penalties.

### 📂 Project Structure
* `Project_2.R`: Complete R source code for data processing, GARCH/Copula fitting, Monte Carlo simulation, and backtesting.
* `Report.pdf`: Detailed presentation slides summarizing the mathematical theory, model selection tables, and visual results.

### 💻 Tech Stack
* **Language**: R
* **Key Libraries**: `rugarch` (GARCH modeling), `rmgarch` (DCC/Multivariate), `copula` (Dependence), `ggplot2` (Visualization), `tseries`, `zoo`.

### 🚀 How to Run
1. Ensure you have the required R packages installed:
   ```r
   install.packages(c("rugarch", "rmgarch", "copula", "ggplot2", "zoo", "xts", "knitr", "moments"))
