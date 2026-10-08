# SSVI / LSSVI learning project

This project constructs an implied-volatility surface from assumed option
quotes using the SSVI parameterization.

Formula:

w(k,theta) = theta/2 * [1 + rho*phi(theta)*k
                         + sqrt((phi(theta)*k + rho)^2 + 1-rho^2)]

phi(theta) = eta / [theta^gamma * (1+theta)^(1-gamma)]

where k=log(K/F) and w=sigma^2*T.

## Files

- data/assumed_option_quotes.csv — synthetic option market data
- calibrate_ssvi.py — nonlinear SSVI calibration
- plot_surface.py — smile and 3-D surface construction
- requirements.txt — dependencies

## Run

pip install -r requirements.txt
python calibrate_ssvi.py
python plot_surface.py

## Important interview point

Newton-Raphson is normally used to invert option prices to implied
volatility. SSVI is then a parametric model for the total-variance surface.
Its parameters are obtained through nonlinear calibration, not by choosing
between Newton-Raphson and cubic splines.

The synthetic data were generated with approximately:
rho=-0.72, eta=0.72, gamma=0.50.

For production data, you would additionally handle bid/ask spreads,
forward/discount construction, vega/liquidity weighting, bad quotes,
butterfly arbitrage, calendar arbitrage, wing extrapolation, and expiry
interpolation.
