Variable Interpretation

Low Rank parameter:
    uvinit: struct.
        uvinit.U: d x R matrix
        uvinit.V: B x R matrix
        uvinit.sigmaU: R x R x d tensor, covariance matrix.
        uvinit.sigmaV: R x R x B tensor, covariance matrix.

MoG Parameter:
    param: struct.
        %%%%%%%%%%%%%%%%%%%%%q(Z)%%%%%%%%%%%%%%%%%%%%%
        param.rau: d x T x B tensor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%q(C)%%%%%%%%%%%%%%%%%%%%%
        param.pphi: T x K x B tensor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%q(\pi^{'})%%%%%%%%%%%%%%%%%%
        param.rrPi1: T x B matrix
        param.rrPi2: T x B matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%q(\beta^{'})%%%%%%%%%%%%%%%%%
        param.ssBeta1: K x 1 vector
        param.ssBeta2: K x 1 vector
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%q(\gamma)%%%%%%%%%%%%%%%%%%%
        param.c:scalar
        param.d:scalar
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%q(\alpha)%%%%%%%%%%%%%%%%%%%
        param.e: B x 1 vector
        param.f: B x 1 vector
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%q(\lambda)%%%%%%%%%%%%%%%%%%
        param.p: 1 scalar
        param.q: R x 1 vector
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    exValue: expectation value of different parameters
        exValue.errSquare: E[(Y-UV')^2]
        %%%%%%%%%%%%%%%%%%q(\lambda)%%%%%%%%%%%%%%%%%%
        exValue.lambda: E[\lambda]
        exValue.logLamda: E[ln(\lambda)]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%q(\alpha)%%%%%%%%%%%%%%%%%%%
        exValue.alpha: E[\alpha], 1 X B vector
        exValue.logAlpha: E[ln(\alpha)], 1 X B vector
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%q(\gamma)%%%%%%%%%%%%%%%%%%%
        exValue.ggamma: E[\gamma]
        exValue.logGamma: E[ln(\gamma)]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%q(\xi)%%%%%%%%%%%%%%%%%%%%%%
        exValue.xi: E[\xi], K x 1 vector
        exValue.logXi: E[log(\xi)], K x 1 vector
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%q(\pi^{'})%%%%%%%%%%%%%%%%%%
        exValue.piPie     : E[\pi^{'}], T x B matrix
        exValue.logPiPie  : E[\ln(\pi^{'})], T x B matrix
        exValue.logPiPie1 : E[\ln(1-\pi^{'})], T x B matrix
        exValue.logPi     : E[\ln(\pi)], T x B matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%q(\beta^{'})%%%%%%%%%%%%%%%%
        exValue.betaPie     : E[\beta^{'}], K x 1 vector
        exValue.logBetaPie1 : E[\ln(\beta^{'})], K x 1 vector
        exValue.logBetaPie  : E[\ln(1-\beta^{'})], K x 1 matrix
        exValue.logBeta     : E[\ln(\beta)], K x 1 matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%q(U)%%%%%%%%%%%%%%%%%%%%
        exValue.UrUr: E[U_{\cdot r}'*E[U_{\cdot r}]]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%q(V)%%%%%%%%%%%%%%%%%%%%
        exValue.UrUr: E[V_{\cdot r}'*E[V_{\cdot r}]]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


