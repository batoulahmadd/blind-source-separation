"""
Embedded Python Blocks:

Each time this file is saved, GRC will instantiate the first class it finds
to get ports and parameters of your block. The arguments to __init__  will
be the parameters. All of them are required to have default values!
"""

import numpy as np
from gnuradio import gr


class blk(gr.sync_block):  # other base classes are basic_block, decim_block, interp_block
    """Embedded Python Block example - a simple multiply const"""

    def __init__(self):  # only default arguments here
        """arguments to this function show up as parameters in GRC"""
        gr.sync_block.__init__(
            self,
            name='correlation Block',   # will show up in GRC
            in_sig=[(np.float32),(np.float32)],  # Two input signals
            out_sig=[]  # No output signal
        )


    def work(self, input_items, output_items):
        signal1 = input_items[0]
        signal2 = input_items[1]
        q=1
        # Calculate the correlation coefficient
        correlation_matrix = np.corrcoef(signal1, signal2)
        correlation_coefficient = correlation_matrix[0, 1]
        print(correlation_coefficient)

        #return len(input_items[0])
