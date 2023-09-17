
import numpy as np
from scipy import linalg


class lst_scalar():
    def __init__(self, N:int, bf:int) -> None:
        """Initialize the LST Analysis

        Args:
            N (int): number of Chebyshev polynomials
            bf (int): set = 1 for Couette, 2 for Poiseuille, 3 for quiescent
            
        Returns:

        """
        # Input parameters 
        self.N       = N # number of Chebyshev polynomials

        R  = 180         # Reynolds number
        kx = 1           # streamwise wavenumber
        kz = 0           # spanwise wavenumber
        Ri = 0.0        # Richardson number
        Pr = 0.71        # Prantl number
        Ra = 8*R*R/Pr*Ri # Rayleigh number
        
        # Set up grid and differentiation matrices
        self.y_phys                = np.cos(np.linspace(0, np.pi, N))[:, np.newaxis]   # Generate Chebyshev grid for base flow solver
        self.D0, self.D1, self.D2, self.D3, self.D4 = self.dmat(N)   # Chebyshev polynomials and derivatives at the Gauss points
        
        
        # Find the base flow
        [U,Up,Upp,T,Tp] = self.bounded_base(self.y_phys,N,bf)
        
        # Find eigenvalues of stability operators
        A, B = self.Operator(kx,kz,R,Pr,Ri,U,Up,Upp,Tp,self.D0, self.D1, self.D2, self.D4)
        
        # find eigenvalues
        # omega, q = linalg.eig(A, B)
        omega, q = np.linalg.eig(np.linalg.inv(B)@A)


        omega = 1j*omega   # eigenvalues omega in vector form

        # remove bad eigenvalues
        sp = np.logical_and(abs(omega)>1e-10, abs(omega)<50)

        self.omega = omega[sp]
        self.q = q[:, sp]
        
        self.Temperature = self.D0 @ self.q
        
        return
        
    def dmat(self, N):
        num = N-1
        D0 = np.cos(np.arange(N)[np.newaxis, :]  * np.pi * np.arange(N)[:, np.newaxis] / num )
        
        # create higher derivative matrices
        D1 = np.concatenate((np.zeros(shape=(N, 1)), D0[:, 0][:, np.newaxis], 4*D0[:, 1][:, np.newaxis]), axis=1)
        D2 = np.concatenate((np.zeros(shape=(N, 2)),                          4*D1[:, 1][:, np.newaxis]), axis=1)
        D3 = np.zeros(shape=(N, 3))
        D4 = np.zeros(shape=(N, 3))



        for j in range(3, N):
            D1= np.concatenate((D1, 2*j*D0[:, j-1][:, np.newaxis]+j*D1[:, j-2][:, np.newaxis]/(j-2)), axis=1)
            D2= np.concatenate((D2, 2*j*D1[:, j-1][:, np.newaxis]+j*D2[:, j-2][:, np.newaxis]/(j-2)), axis=1)
            D3= np.concatenate((D3, 2*j*D2[:, j-1][:, np.newaxis]+j*D3[:, j-2][:, np.newaxis]/(j-2)), axis=1)
            D4= np.concatenate((D4, 2*j*D3[:, j-1][:, np.newaxis]+j*D4[:, j-2][:, np.newaxis]/(j-2)), axis=1)
            
        return D0, D1, D2, D3, D4
    
    def bounded_base(self, y_phys,N,bf):
        if bf == 1: # Couette flow
    
            U   = y_phys
            Up  = np.ones(shape=(N, 1))
            Upp = np.zeros(shape=(N, 1))
            T   = y_phys
            Tp  = np.zeros(shape=(N, 1))
            
        elif bf == 2 : # Poiseuille flow
            
            U   = (1 - y_phys**2)
            Up  = -2*y_phys
            Upp = -2*np.ones(shape=(N, 1))
            T   = np.zeros(shape=(N, 1))
            Tp  = np.zeros(shape=(N, 1))
            
        elif bf == 3: # quiescent flow
            
            U   = np.zeros(shape=(N, 1))
            Up  = np.zeros(shape=(N, 1))
            Upp = np.zeros(shape=(N, 1))
            T   = y_phys
            Tp  = np.ones(shape=(N, 1))
        else:
            print('Need to select bf = 1 or 2')
            
            
        return U,Up,Upp,T,Tp
    
    def Operator(self, kx,kz,R,Pr,Ri,U,Up,Upp,Tp,D0,D1,D2,D4):
        N = self.N
        
        k2 = kx**2 + kz**2 
        M  =  np.ones(shape=(1, N))
        er = -200*1j 

        A = -1j*kx*(U@M)*D0 + (1/R/Pr)*(D2-k2*D0)
        B = D0


        A[0, :] = 0
        A[N-1, :] = 0

        A[0, :] = er*D0[0, :]
        A[N-1, :] = er*D0[N-1, :]
        
        return A, B
    
    def visual(self, axis, eigvals, eigvecs, marker_types = [dict(size=12, symbol="circle", ), dict(size=8, symbol="x", )]):
        from dash import Dash, html, dcc, Input, Output, callback
        import pandas as pd
        import plotly.express as px
        import plotly.graph_objects as go
        
        # Prepare LST result in df1
        df1 = pd.DataFrame(
        dict(
            # eigvals = channel_lst.omega,
            eigvals_real = np.real(self.omega),
            eigvals_imag = np.imag(self.omega),
            eigvecs = [self.Temperature[:, ind] for ind in range(self.omega.shape[0])],
        )
        )
        y_phys = self.y_phys[:, 0]

        df1["group"] = "LST"
        
        
        # Prepare solver result in df2
        Nz, Npoint = eigvecs.shape
        df2 = pd.DataFrame(
        dict(
            eigvals_real = np.real(eigvals).reshape((Npoint)),
            eigvals_imag = np.imag(eigvals).reshape((Npoint)),
            eigvecs = [eigvecs[:, ind] for ind in range(Npoint)],
        )
        )
        df2["group"] = "solver"
        
        y_set = [y_phys, axis]
        
        # Assemble
        frames = [df1, df2]
        df = pd.concat(frames, ignore_index=True)
        external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

        app = Dash(__name__, external_stylesheets=external_stylesheets)


        app.layout = html.Div([
            html.Div([

                html.Div([
                    "Layer Number: ",
                    dcc.Input(id='layer_number', value=384, type='number')
                ], style={'width': '49%', 'display': 'inline-block', 'padding': '0 20'}),
            ]),
            html.Br(),

            html.Div([
                html.Div([
                    html.Label('Display result'),
                    dcc.Checklist(id="display_classes", options=df["group"].unique(),
                                value = []
                    ),
                ], style={'width': '49%', 'display': 'inline-block', 'padding': '0 20'}),

                
                html.Div(children=[
                    html.Div(children=[
                        "LST Index: ",
                        dcc.Input(id='lst_index', value=0, type='number'),
                        html.Div(id="group1_eigenvalues" )
                    ], ),
                    
                    html.Div([
                        "solver Index: ",
                        dcc.Input(id='result_index', value=0, type='number'),
                        html.Div(id="group2_eigenvalues" )
                    ],),
                ], style={'width': '24%', 'display': 'inline-block'}),
            ]),
            
            html.Div([

                html.Div([
                    dcc.Graph(
                        id='eigenvals_plot',
                        clickData={'points': [{'hovertext': 0, 'curveNumber' : 0}]}
                    )
                ], style={'width': '49%', 'display': 'inline-block', 'padding': '0 20'}),


                html.Div([
                    dcc.Graph(id='eigenvector_plot',)
                ], style={'display': 'inline-block', 'width': '49%'})
            ]),
        ])

        @callback(
            Output('eigenvals_plot', 'figure'),
            Input('layer_number', 'value'),
            Input("display_classes", "value")
        )
        def update_graph(layer_number, display_class):
            # fig = px.scatter(df, x="eigvals_real", y = "eigvals_imag", color="group", symbol="group", hover_name=df.index)
            
            fig = go.Figure()

            # Add traces
            for ind, cls in enumerate(df["group"].unique()):
                if cls in display_class:
                    fig.add_trace(go.Scatter(x=df[df["group"]==cls]["eigvals_real"], y = df[df["group"]==cls]["eigvals_imag"], hovertext=df[df["group"]==cls].index, opacity=0.5, marker=marker_types[ind],
                                        mode='markers',
                                        name=df["group"].unique()[ind]))
            
            fig.update_layout(template="simple_white", height=450, margin={'l': 20, 'b': 30, 'r': 10, 't': 10})
                    
            return fig


        @callback(
            Output('eigenvector_plot', 'figure'),
            Output('lst_index', 'value'),
            Output('result_index', 'value'),
            Output("group1_eigenvalues", "children"),
            Output("group2_eigenvalues", "children"),
            Input('eigenvals_plot', 'clickData'),
            Input('layer_number', 'value'),
            Input('lst_index', 'value'),
            Input('result_index', 'value'),
        )
        def update_eigenmode(clickData, layer_number, lst_ind, result_ind):
            norm_eigvec = lambda eigvecs: np.abs(eigvecs)/np.max(np.abs(eigvecs))
            
            
            print(clickData)
            eig_ind = clickData["points"][0]['hovertext']
            group = clickData["points"][0]["curveNumber"] # 0 - LST, 1 - result
            
            ind_set = [lst_ind, result_ind]
            if eig_ind in df[df["group"]==df["group"].unique()[0]].index:
                ind_set[0] = eig_ind
            if eig_ind in df[df["group"]==df["group"].unique()[1]].index:
                ind_set[1] = eig_ind


            # fig = px.line(x= y_phys, y= np.abs(df['eigvecs'][lst_ind]))
            mode_types = ["markers", "lines"]
            
            fig = go.Figure()

            # Add traces
            for ind, cls in enumerate(df["group"].unique()):
                    fig.add_trace(go.Scatter(x=y_set[ind], y=norm_eigvec(np.abs(df['eigvecs'][ind_set[ind]])), marker=marker_types[ind],
                                        mode=mode_types[ind],
                                        name=df["group"].unique()[ind]))
                # fig.add_trace(go.Scatter(x=result_y, y=np.abs(df['eigvecs'][result_ind]),
                #                     mode='lines',
                #                     name=df["group"].unique()[1]))
                
            fig.update_layout(template="simple_white", height=450, margin={'l': 20, 'b': 30, 'r': 10, 't': 10})
            
            group1_eigvals = "%f + %f i" %(df['eigvals_real'][ind_set[0]], df['eigvals_imag'][ind_set[0]])
            group2_eigvals = "%f + %f i" %(df['eigvals_real'][ind_set[1]], df['eigvals_imag'][ind_set[1]])



            return fig, ind_set[0], ind_set[1], group1_eigvals, group2_eigvals
        return app