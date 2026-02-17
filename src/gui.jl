# ============================================================================
# GUI — Web-based Dashboard for Phase Transition Simulation (WGLMakie + Bonito)
# ============================================================================

using WGLMakie
using WGLMakie: Bonito
using Observables
using Markdown

"""
Launch the web-based interactive Ising model phase transition simulation GUI.
"""
function launch_gui(; port=8080)
    # ── Simulation State ────────────────────────────────────────────────
    lattice_type = Observable("Square")
    lattice_size = Observable(32)
    temperature  = Observable(2.25)
    coupling_J   = Observable(1.0)
    field_h      = Observable(0.0)
    bc_type      = Observable("Periodic")
    algo_type    = Observable("Metropolis")
    running      = Observable(false)
    sweep_mode   = Observable(false)

    # Current lattice & model
    current_lat   = Observable{AbstractLattice}(SquareLattice(32))
    current_model = Observable(IsingModel(1.0, 0.0))

    # Live data (Observables for plotting)
    spin_data        = Observable(zeros(Float32, 32, 32))
    energy_history   = Observable(Float64[])
    mag_history      = Observable(Float64[])
    step_indices     = Observable(Int[])
    
    # Sweep results
    sweep_T     = Observable(Float64[])
    sweep_E     = Observable(Float64[])
    sweep_absM  = Observable(Float64[])
    sweep_Cv    = Observable(Float64[])
    sweep_chi   = Observable(Float64[])
    sweep_U4    = Observable(Float64[])

    # Analysis data
    autocorr_data   = Observable(Float64[])
    autocorr_lags   = Observable(Float64[])
    corr_func_data  = Observable(Float64[])
    corr_func_r     = Observable(Float64[])
    
    # Derived scalar values for display
    tau_int_val     = Observable(0.0)
    corr_length_val = Observable(0.0)
    status_text     = Observable("Ready")

    # Sweep parameters
    T_start = Observable(1.0)
    T_end   = Observable(4.0)

    # ── Helper Functions ────────────────────────────────────────────────

    function rebuild_lattice!()
        bc = bc_type[] == "Periodic" ? Periodic : Open
        lat = create_lattice(lattice_type[], lattice_size[]; bc=bc)
        current_lat[] = lat
        current_model[] = IsingModel(coupling_J[], field_h[])
        
        # Reset histories
        energy_history[] = Float64[]
        mag_history[] = Float64[]
        step_indices[] = Int[]
        
        spin_data[] = Float32.(spin_matrix(lat))
    end

    function update_analysis!()
        e_hist = energy_history[]
        if length(e_hist) > 20
            ac = autocorrelation(e_hist; max_lag=min(200, length(e_hist) ÷ 4))
            autocorr_data[] = ac
            autocorr_lags[] = collect(Float64, 0:length(ac)-1)
            tau_int_val[] = integrated_autocorrelation_time(e_hist)
        end

        lat = current_lat[]
        G = spatial_correlation(lat)
        corr_func_data[] = G
        corr_func_r[] = collect(Float64, 0:length(G)-1)
        corr_length_val[] = correlation_length(lat)
    end

    # ── Bonito App ──────────────────────────────────────────────────────

    app = Bonito.App() do session
        # Use a dark theme for plots
        set_theme!(theme_dark())

        # -- Controls --
        
        # Lattices
        lat_dropdown = Bonito.Dropdown(LATTICE_NAMES; index=1)
        on(v -> lattice_type[] = v, lat_dropdown.value)
        
        size_slider = Bonito.Slider(8:4:64; value=32)
        on(v -> lattice_size[] = Int(v), size_slider.value)
        
        bc_dropdown = Bonito.Dropdown(["Periodic", "Open"]; index=1)
        on(v -> bc_type[] = v, bc_dropdown.value)

        # Physics
        temp_slider = Bonito.Slider(0.1:0.05:6.0; value=2.25)
        on(v -> temperature[] = v, temp_slider.value)
        
        j_slider = Bonito.Slider(-3.0:0.1:3.0; value=1.0)
        on(v -> (coupling_J[] = v; current_model[] = IsingModel(v, field_h[])), j_slider.value)
        
        h_slider = Bonito.Slider(-2.0:0.1:2.0; value=0.0)
        on(v -> (field_h[] = v; current_model[] = IsingModel(coupling_J[], v)), h_slider.value)

        algo_dropdown = Bonito.Dropdown(["Metropolis", "Wolff"]; index=1)
        on(v -> algo_type[] = v, algo_dropdown.value)

        # Buttons
        start_btn = Bonito.Button("Start"; style=Bonito.Styles(Bonito.CSS("background-color" => "#4CAF50"), Bonito.CSS("color" => "white")))
        stop_btn  = Bonito.Button("Stop"; style=Bonito.Styles(Bonito.CSS("background-color" => "#f44336"), Bonito.CSS("color" => "white")))
        reset_btn = Bonito.Button("Reset")
        sweep_btn = Bonito.Button("Temperature Sweep"; style=Bonito.Styles(Bonito.CSS("background-color" => "#2196F3"), Bonito.CSS("color" => "white")))

        # Button Logic
        on(start_btn.value) do _
            if !running[]
                rebuild_lattice!()
                running[] = true
                sweep_mode[] = false
                status_text[] = "Running..."
                
                @async begin
                    try
                        lat = current_lat[]
                        mod = current_model[]
                        algo = algo_type[] == "Wolff" ? :wolff : :metropolis
                        step = 0
                        while running[]
                            β = 1.0 / temperature[]
                            mc_sweep!(lat, mod, β, algo)
                            step += 1
                            if step % 5 == 0
                                e = energy_per_spin(mod, lat)
                                m = magnetization_per_spin(lat)
                                push!(energy_history[], e)
                                push!(mag_history[], m)
                                push!(step_indices[], step)
                                
                                if length(energy_history[]) > 1000
                                    deleteat!(energy_history[], 1:100)
                                    deleteat!(mag_history[], 1:100)
                                    deleteat!(step_indices[], 1:100)
                                end
                                
                                notify(energy_history); notify(mag_history); notify(step_indices)
                                spin_data[] = Float32.(spin_matrix(lat))
                                
                                if step % 50 == 0
                                    update_analysis!()
                                end
                                sleep(0.001)
                            end
                        end
                    catch e
                        println("Error: $e")
                        running[] = false
                    end
                end
            end
        end

        on(stop_btn.value) do _
            running[] = false
            status_text[] = "Stopped"
            update_analysis!()
        end

        on(reset_btn.value) do _
            running[] = false
            rebuild_lattice!()
            status_text[] = "Reset"
        end

        on(sweep_btn.value) do _
            if !running[]
                running[] = true
                sweep_mode[] = true
                status_text[] = "Sweeping..."
                
                @async begin
                    try
                        bc = bc_type[] == "Periodic" ? Periodic : Open
                        algo = algo_type[] == "Wolff" ? :wolff : :metropolis
                        mod = IsingModel(coupling_J[], field_h[])
                        L = lattice_size[]
                        lat_name = lattice_type[]
                        
                        trange = range(T_start[], T_end[], length=30)
                        
                        Ts, Es, aMs, Cvs, chis, U4s = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]
                        
                        for T in trange
                            if !running[] break end
                            β = 1.0 / T
                            lat = create_lattice(lat_name, L; bc=bc)
                            for _ in 1:1000; mc_sweep!(lat, mod, β, algo); end # therm
                            
                            e_sum=0.0; e2_sum=0.0; m2_sum=0.0; m4_sum=0.0; am_sum=0.0
                            N = nsites(lat)
                            measure_steps = 1000
                            
                            for _ in 1:measure_steps
                                mc_sweep!(lat, mod, β, algo)
                                e = energy(mod, lat)/N; m = magnetization(lat)/N
                                e_sum += e; e2_sum += e^2
                                m2_sum += m^2; m4_sum += m^4; am_sum += abs(m)
                            end
                            
                            push!(Ts, T)
                            push!(Es, e_sum/measure_steps)
                            push!(aMs, am_sum/measure_steps)
                            
                            e_avg = e_sum/measure_steps
                            e2_avg = e2_sum/measure_steps
                            cv = β^2 * N * (e2_avg - e_avg^2)
                            push!(Cvs, cv)
                            
                            m2_avg = m2_sum/measure_steps
                            am_avg = am_sum/measure_steps
                            chi = β * N * (m2_avg - am_avg^2)
                             push!(chis, chi)

                            m4_avg = m4_sum / measure_steps
                            u4 = m2_avg ≈ 0.0 ? 0.0 : 1.0 - m4_avg / (3.0 * m2_avg^2)
                            push!(U4s, u4)
                            
                            sweep_T[] = copy(Ts); sweep_E[] = copy(Es); sweep_absM[] = copy(aMs)
                            sweep_Cv[] = copy(Cvs); sweep_chi[] = copy(chis); sweep_U4[] = copy(U4s)
                            
                            spin_data[] = Float32.(spin_matrix(lat))
                            status_text[] = "Sweep T=$T"
                            sleep(0.01)
                        end
                        running[] = false
                        status_text[] = "Sweep Done"
                    catch e
                         println("Error in sweep: $e")
                         running[] = false
                    end
                end
            end
        end

        # -- GUI Layout --
        
        # Styles
        card_style = "background-color: #1e1e24; padding: 15px; border-radius: 8px; margin: 5px; color: white;"
        
        # Plots
        fig = Figure(size=(1000, 800), backgroundcolor="#1e1e24")
        
        ax_spin = Axis(fig[1, 1:2], title="Spin Configuration")
        heatmap!(ax_spin, spin_data, colormap=:RdBu, colorrange=(-1, 1))
        
        ax_en = Axis(fig[2, 1], title="Energy", xlabel="Step")
        lines!(ax_en, step_indices, energy_history, color=:cyan)
        
        ax_mag = Axis(fig[2, 2], title="Magnetization", xlabel="Step")
        lines!(ax_mag, step_indices, mag_history, color=:orange)
        
        ax_cv = Axis(fig[3, 1], title="Specific Heat vs T")
        scatterlines!(ax_cv, sweep_T, sweep_Cv, color=:red)
        
        ax_chi = Axis(fig[3, 2], title="Susceptibility vs T")
        scatterlines!(ax_chi, sweep_T, sweep_chi, color=:green)
        
        # Layout construction
        layout = Bonito.DOM.div(
            Bonito.DOM.h1("Ising Model Simulator", style="color: white; text-align: center;"),
            Bonito.DOM.div(
                Bonito.DOM.div(
                    Bonito.DOM.h3("Controls"),
                    Bonito.DOM.label("Lattice:"), lat_dropdown,
                    Bonito.DOM.label("Size:"), size_slider,
                    Bonito.DOM.label("Boundary:"), bc_dropdown,
                    Bonito.DOM.hr(),
                    Bonito.DOM.label("Temperature:"), temp_slider,
                    Bonito.DOM.label("Coupling J:"), j_slider,
                    Bonito.DOM.label("Field h:"), h_slider,
                    Bonito.DOM.label("Algorithm:"), algo_dropdown,
                    Bonito.DOM.hr(),
                    Bonito.DOM.div(start_btn, stop_btn, reset_btn, sweep_btn, style="display: flex; gap: 5px; flex-wrap: wrap;"),
                    Bonito.DOM.p(status_text, style="color: yellow; font-weight: bold; margin-top: 10px;"),
                    style=card_style * "width: 250px; flex-shrink: 0;"
                ),
                Bonito.DOM.div(
                    fig,
                    style=card_style * "flex-grow: 1;"
                ),
                style="display: flex; flex-direction: row; gap: 10px;"
            ),
            style="font-family: sans-serif; background-color: #121212; min-height: 100vh; padding: 20px;"
        )

        return layout
    end

    # Launch
    Bonito.Server(app, "0.0.0.0", port)
    println("Server running at http://localhost:$port")
    wait(Condition()) # Keep alive
end
