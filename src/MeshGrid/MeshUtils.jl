## Node position for each element
mutable struct NodalCoordinatesInElement{T<:Array}
    x::T
    y::T
    z::T
end
function NodePosition3D(mesh::AbstractMesh)
    table = size(mesh.IEN)
    EN_x = zeros(Float64, table) #počet ele = počet Gausspointů
    EN_y = zeros(Float64, table)
    EN_z = zeros(Float64, table)

    for i = 1:(length(mesh.IEN[:, 1]))
        for j = 1:(length(mesh.IEN[1, :]))
            ID_n = Int(mesh.IEN[i, j])
            x = mesh.X[1, ID_n]
            y = mesh.X[2, ID_n]
            z = mesh.X[3, ID_n]
            EN_x[i, j] = x
            EN_y[i, j] = y
            EN_z[i, j] = z
        end
    end
    EN = NodalCoordinatesInElement(EN_x, EN_y, EN_z)

    return EN #Výstupem jsou tři pole - souřadnice uzlů daného elementu (pozice řádku = ID elemetu)
end


# Geometric centres of elements (GPs for first order elements)
function GeometricCentre(mesh::AbstractMesh, EN::NodalCoordinatesInElement)
    Centre = zeros(mesh.nel, length(mesh.X[:, 1]))
    for i = 1:mesh.nel
        G_x = mean(EN.x[:, i]) # G_x = dot(EN_x[:, i], N)
        G_y = mean(EN.y[:, i]) # G_y = dot(EN_y[:, i], N)
        G_z = mean(EN.z[:, i]) # G_z = dot(EN_z[:, i], N)
        Centre[i, :] = [G_x, G_y, G_z]
    end
    return Centre
end
