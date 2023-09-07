using MAT, LinearAlgebra
using Dash, PlotlyJS, DataFrames

N = 384

lst_vars = matread("/Users/user/Documents/Projects/python_util/test/data/scalar_LST.mat")
scalar_eigvals = eigvals(inv(lst_vars["B"]) * lst_vars["A"] ) * 1im
scalar_eigvecs = eigvecs(inv(lst_vars["B"]) * lst_vars["A"] )


# remove bad eigenvalues
scalar_sp = (abs.(scalar_eigvals).>1e-10) .*  (abs.(scalar_eigvals).<50)

scalar_eigvals = scalar_eigvals[scalar_sp]
scalar_eigvecs = lst_vars["D0p"] * scalar_eigvecs[:, scalar_sp]


df = DataFrame(
    eigvals = scalar_eigvals,
    eigvecs = [scalar_eigvecs[:, i] for i in range(1, size(scalar_eigvecs)[2])],
    eig_ind = range(1, length(scalar_eigvals))
)



app = dash()

app.layout = html_div() do
    html_div() do
        html_div(
            children = [
                "Layer Number: ",
                dcc_input(
                    id="layer_number", 
                    value=384, 
                    type="number"
                ),
            ], 
            style=(width = "48%", display = "inline-block", padding = "0 20"),
        ),
        
        
        html_div(
            children = [
                "Point Index: ",
                dcc_input(id="point_index", value=0, type="number"),
            ], 
            style=(width = "48%", display = "inline-block"),
        )
    end,

    html_div() do
        html_div(
            children = [
                dcc_graph(
                    id="eigenvals_plot",
                    clickData=Dict("points"=> [Dict("pointNumber"=> 1)])
                ),
            ], 
            style=(width = "48%", display = "inline-block", padding = "0 20"),
        ),
            

        html_div(
            children = [
                dcc_graph(
                    id="eigenvector_plot",            
                ),
            ], 
            style=(width = "48%", display = "inline-block"),
        )
    end
end

callback!(
    app,
    Output("eigenvals_plot", "figure"),
    Input("layer_number", "value")
) do layer_number
    fig = plot(scatter(x=real(df.eigvals), y = imag(df.eigvals), mode="markers"))
    
    return fig
end

callback!(
    app,
    Output("eigenvector_plot", "figure"),
    Output("point_index", "value"),
    Input("eigenvals_plot", "clickData"),
    Input("layer_number", "value")
) do clickData, layer_number
    y_phys = cos.(LinRange(0, pi, layer_number))
    eig_ind = clickData["points"][1]["pointNumber"]
    print(eig_ind)
    fig = plot(y_phys, abs.(df.eigvecs[eig_ind]))
    return fig, eig_ind
end

run_server(app, "127.0.0.1", 8088, debug=true)