"""
    export_structured_field_to_vti(
        filename::String,
        grid::Grid,
        values::Vector{<:AbstractFloat},
        value_label::String,
        smooth::Union{Int,Nothing} = nothing)

Export skalárního pole definovaného na strukturované mřížce do VTI souboru (VTK ImageData).
VTI formát je optimalizovaný pro pravidelné mřížky s konstantním krokem ve všech směrech.

# Argumenty
- `filename::String`: Název výstupního VTI souboru (bez přípony, automaticky se přidá .vti)
- `grid::Grid`: Strukturovaná mřížka, na které jsou definovány hodnoty
- `values::Vector{<:AbstractFloat}`: Skalární hodnoty v každém bodě mřížky
- `value_label::String`: Popisek skalárního pole (např. "distance")
- `smooth::Union{Int,Nothing} = nothing`: Volitelný faktor vyhlazení pro generování jemnější mřížky

# Návratová hodnota
- Odkaz na vytvořený VTK objekt
"""
function exportSdfToVTI(
    filename::String,
    grid::Grid,
    values::Vector{<:AbstractFloat},
    value_label::String,
    smooth::Union{Int,Nothing} = nothing
)
    # Výpočet dimenzí, počátku a kroku mřížky pro VTI formát
    dimensions = isnothing(smooth) ? grid.N .+ 1 : (grid.N * smooth) .+ 1
    origin = Float64.(grid.AABB_min)
    
    # Pokud je spacing skalár, potřebujeme ho převést na vektor pro každou dimenzi
    if isa(grid.cell_size, Number)
        spacing = isnothing(smooth) ? 
                 [grid.cell_size, grid.cell_size, grid.cell_size] :
                 [grid.cell_size/smooth, grid.cell_size/smooth, grid.cell_size/smooth]
    else
        spacing = isnothing(smooth) ? 
                 grid.cell_size : 
                 grid.cell_size ./ smooth
    end
    
    # Ověření správné délky vektoru hodnot
    if length(values) != prod(dimensions)
        error("Délka vektoru hodnot ($(length(values))) neodpovídá dimenzím mřížky ($(prod(dimensions))).")
    end
    
    # Generování souřadnic podél každé osy pro vytvoření strukturované mřížky
    x = range(origin[1], length=dimensions[1], step=spacing[1])
    y = range(origin[2], length=dimensions[2], step=spacing[2])
    z = range(origin[3], length=dimensions[3], step=spacing[3])
    
    # Vytvoření VTI strukturované mřížky - použití správné signatury
    vtk_file = vtk_grid(filename, x, y, z)
    
    # Přetvarování hodnot do 3D pole pro WriteVTK
    values_3d = reshape(values, dimensions...)
    
    # Přidání dat skalárního pole
    vtk_point_data(vtk_file, values_3d, value_label)
    
    # Uložení VTI souboru
    vtk_save(vtk_file)
    
    return vtk_file
end
