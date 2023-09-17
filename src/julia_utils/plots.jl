using Plots

# Define a callback function to capture click data
function click_callback(event)
    if event == Plots.MouseEvent(:click, 1)
        println("Left-clicked at position: ", current_point())
    elseif event == Plots.MouseEvent(:click, 3)
        println("Right-clicked at position: ", current_point())
    end
end

# Create a simple scatter plot
scatter([1, 2, 3, 4, 5], [5, 4, 3, 2, 1], label="Points")

# Attach the click callback to the plot
on(1, click_callback)

# Display the plot
display(1)
