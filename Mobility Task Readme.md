This is a great task, and it's a significant step up from the previous simulation. In Task 1/2, we had a simple **AWGN** channel, which just adds noise. In this task, we are simulating a **TDL (Tapped Delay Line)** channel, which is a realistic 5G model for a *fading* channel.

A fading channel doesn't just add noise; it actively distorts the signal. It creates echoes (delay spread) and, when you move, it rapidly changes the signal's phase and amplitude (Doppler shift).

Because the channel is so much more complex, our simulation has to be more complex too. We can't just send bits. We have to build a full 5G-style transmitter and receiver.

Here is a full, one-by-one explanation of the entire script.

---

### **Section 1: Master Simulation Parameters**

This section just sets up all the variables for our experiment.

* `K = 54`, `E = 124`, etc.: These are our basic Polar code settings, carried over from the last script.
* `EbNo_fixed_dB = 4.0`: This is a key decision. In Task 1, our X-axis was `EbNo`. In this task, our X-axis is `speed_kmh_vec`. To draw a 2D graph, we *must* fix the signal quality. We chose 4.0 dB because our Task 2 results likely showed this is a challenging "cell-edge" region where the choice of `L` *matters*.
* `L_vec = [2, 4, 8, 16, 32]`: These are the "various options" your professor mentioned. Each `L` value will get its own line on the final graph.
* `speed_kmh_vec = [3, 30, 60, 120]`: This is our new **X-axis**, the "changed setting". It covers everything from walking (3 km/h) to highway speed (120 km/h).
* `maxNumFrames`, `minFrameErrors`: Same as before, just simulation controls for accuracy.
* `carrierFrequency_Hz = 4e9`: This is **critical** for a mobility simulation. The *Doppler shift* (the distortion from movement) is directly proportional to the carrier frequency. 4 GHz is a typical 5G "mid-band" frequency.
* `channelProfile = 'TDL-C'`: This tells the simulator what *kind* of fading to use. 'TDL-C' is a standard 3GPP model for a typical urban, non-line-of-sight (NLOS) environment.

---

### **Section 2: 5G NR PDSCH & TDL Channel Configuration**

This is the biggest new part. Because we have a fading channel, we need to build a realistic 5G signal structure that knows how to *fight* the fading.

* `pdsch = nrPDSCHConfig`: We create a **PDSCH (Physical Downlink Shared Channel)** configuration object. This is the 5G data channel. This object is a blueprint that defines...
    * `NumLayers = 1`: How many antennas (we'll keep it simple).
    * `PRBSet = 0:10`: *Where* in the 5G frequency band we will put our data.
    * `SymbolAllocation = [0 14]`: We'll use the whole 5G "slot" in time.
* `pdsch.DMRS...`: This is the **most important part for a fading channel**. DMRS stands for **Demodulation Reference Signals**, also known as **pilots**.
    * **Analogy:** Imagine trying to read a sign in a funhouse mirror (the fading channel). It's impossible. But what if you tape a perfect, square grid onto the mirror first? By seeing how the grid is distorted, you can *learn* the shape of the mirror and *mentally reverse* the distortion, allowing you to read the sign.
    * The DMRS pilots are that "perfect grid." They are known signals we "sprinkle" into our transmission so the receiver can measure the channel's distortion. 
* `ofdmInfo = nrOFDMInfo(...)`: This defines the fundamental structure of the 5G signal (e.g., its sample rate, which is derived from the channel's 15kHz subcarrier spacing).
* `tdl = nrTDLChannel`: This creates our "enemy." This is the MATLAB object that *simulates* the TDL-C fading channel. We set its sample rate to match our 5G signal.

---

### **Section 3: Simulation Processing Loop (Task 5)**

This is where the actual experiment happens. We have two outer loops that set up the plot, and one inner `while` loop that runs the simulation for each point.

**Outer Loops (`i_L` and `i_speed`):**
These are simple. The first loop picks a list size (`L=2`), and the second loop iterates through each speed (`3`, `30`, `60`, `120`).

**Inside the Speed Loop:**
* `speed_mps = speed_kmh / 3.6`: Convert km/h to m/s.
* `dopplerShift_Hz = ...`: This is the core physics. This formula converts the UE's speed (in m/s) and the carrier frequency (in Hz) into the **Maximum Doppler Shift** (in Hz). This value tells the channel simulator *how fast* to change the fading. Higher speed = higher Doppler = faster, more difficult channel.
* `tdl.MaximumDopplerShift = ...`: We *configure* our channel object (`tdl`) with this new Doppler shift for the current simulation point.
* `[~, pdschInfo] = nrPDSCHIndices(pdsch)`: A helper function. We ask the toolbox, "Given our PDSCH config, how many coded bits (`E_actual`) *actually fit* on the resource grid?" This is more accurate than our old, fixed `E=124`.
* `codeRate = K/E_actual`: We calculate the *true* code rate.
* `noiseVar = ...`: Same as before. We calculate the noise power needed for our fixed `EbNo_fixed_dB` value.
* `reset(tdl)`: We reset the channel's internal state for a fair test at each new speed.

**The `while` Loop (One Frame Transmission):**
This is the step-by-step process of sending one block of data and checking if it failed.

#### **TRANSMITTER (TX)**
* `msg`, `msg_crc`, `rateMatchedBits`: Same as before. We create a message, add a CRC check, and Polar encode it.
* `pdschGrid = nrPDSCH(pdsch, rateMatchedBits)`: This takes our encoded bits and *maps* them onto the correct locations of a blank 2D resource grid.
* `dmrsSymbols = nrPDSCHDMRS(pdsch)`: Create the *pilot* symbols.
* `pdschGrid(dmrsIndices) = dmrsSymbols`: *Pokes* the pilots into the grid. The grid now contains both our data and the known reference signals.
* `txWaveform = nrOFDMModulate(...)`: This converts the 2D grid (frequency/time) into a 1D time-domain *waveform* (a complex radio signal) that is ready to be "sent."

#### **CHANNEL (CH)**
* `[rxWaveform, pathGains] = tdl(txWaveform)`: The waveform is passed through our TDL channel simulator. The output `rxWaveform` is the *distorted, faded* version of our signal.
* `noise = ...`, `rxWaveform_noisy = ...`: We add the AWGN *after* the fading, just like in the real world.

#### **RECEIVER (RX)**
This is the reverse of the TX, but with two new, critical steps.
* `rxGrid = nrOFDMDemodulate(...)`: Converts the 1D distorted waveform *back* into a 2D resource grid. This grid is a mess—it's been distorted by fading and noise.
* `[estChannel, estNoiseVar] = nrChannelEstimate(...)`: **This is Step 1 of the "antidote."** This is the "practical" estimation from your evaluation plan. We tell this function *where* the pilots were (`dmrsIndices`) and *what* they should have been (`dmrsSymbols`). It compares them to what it *actually* received in those locations and creates a "map" (`estChannel`) of the distortion across the *entire* grid. This is our *guess* of what the funhouse mirror looks like.
* `[pdschSymbols, ...], [pdschEstH, ...]`: We *extract* just the data symbols (and their corresponding distortion estimates) from the grid.
* `eqSymbols = nrEqualizeMMSE(...)`: **This is Step 2 of the "antidote."** This is the **MMSE Equalizer**. It takes the distorted data symbols (`pdschSymbols`) and our channel map (`pdschEstH`) and *reverses the distortion*. The output `eqSymbols` is our *best guess* of what the symbols looked like *before* the fading channel.
* `rxLLRs = nrSymbolDemodulate(...)`: Same as before. We turn the "clean" (equalized) symbols into LLRs (soft bits).
* `decodedBits = nrPolarDecode(...)`: Same as before. The Polar decoder (using the current loop's `L` value) tries to recover the message from the LLRs.
* `if any(...)`: We check if the decoded message matches the one we sent. If not, we increment `errorCount`.

---

### **Section 4: Plot Results (Task 5)**

This section just plots the `bler_results` matrix we spent all that time filling.

* `figure; semilogy(...)`: Create a new figure with a logarithmic Y-axis.
* `xlabel('UE Speed (km/h)')`: Our **X-axis** is the UE speed.
* `ylabel('Block Error Rate (BLER)')`: Our **Y-axis** is the BLER.
* `legend(...)`: This adds the legend, where each line is labeled with its List Size (e.g., "L = 2", "L = 4").

**What You Will See:**
The resulting graph will show all the lines "going up" (BLER getting worse) as speed increases. This is because the Doppler shift is making the channel change too fast for the receiver to accurately estimate and correct it.

Your key analysis, as you defined in your goal, will be to look at the *gap between the lines*. You will likely see that the performance gain from `L=32` over `L=8` is *smaller* at 120 km/h than it is at 3 km/h. This would be a great conclusion: "At high speeds, the channel quality degrades so much that the limiting factor is the channel estimation, and the benefit of a more complex decoder is diminished."