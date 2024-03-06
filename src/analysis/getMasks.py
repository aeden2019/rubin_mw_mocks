"""
Functions that return masks of MSTO, BHB, and K-Giant stellar populations.

The populations are defined by Section 3.4 of Lowing 2018: https://arxiv.org/pdf/1406.2320.pdf

"""



def getMSTO(g_minus_r, g_mag):
    """
    To select main sequence turn off stars we use $ 0.2 < g−r < 0.4$ and $M_g > 4$.
    
    Parameters:
    -----------
    g_minus_r: array of the g-r color band
    g_mag: array of absolute magnitudes in the g-band
    
    
    Returns:
    --------
    msto_mask: boolean mask of MSTO stars
    
    """
    
    msto_mask = ((g_minus_r > 0.2) & 
                (g_minus_r < 0.4) & 
                (g_mag > 4))
    
    return msto_mask



def getBHB(g_minus_r, u_minus_g, g_mag):
    """
    To select BHB stars  we use $ 0.98 < u − g < 1.28$, $−0.27 < g − r < −0.06$ and excluding the region $([u−g −0.98]/0.215)2 +
    ([g −r+ 0.06]/0.17)2 < 1$.
    
    Parameters:
    -----------
    g_minus_r: array of the g-r color band
    u_minus_g: array of the u-g color band
    g_mag: array of absolute magnitudes in the g-band
    
    Returns:
    --------
    bhb_mask: boolean mask of BHB stars
    
    """
    
    bhb_mask = ((u_minus_g > 0.98) & 
                (u_minus_g < 1.28) &
                (g_minus_r > -0.27) & 
                (g_minus_r < -0.06) &
                ((((u_minus_g - 0.215) / 0.215) * 2 + 
                  ((g_minus_r + 0.06) / 0.17) * 2) >= 1))
    
    return bhb_mask



def getKGiant(g_minus_r, u_minus_g, g_mag, feh):
    """
    In order to isolate K giants in our stellar halo models we have used the colour cuts $0.5 < g − r < 1.3$ and $0.5 <
    u−g < 3.5$ from Xue et al. (2014), along with their proposed empirical polynomial relation between (g −r) and [Fe/H] to remove
    red horizontal branch and red-clump (RC) giants. All stars with $(g − r) > 0.087 [Fe/H]2 + 0.39 [Fe/H] + 0.96$ is excluded
    from the selection. An additional cut of $M_g < 4$ removes faint dwarf stars of the same colours
    
    Parameters:
    -----------
    g_minus_r: array of the g-r color band
    u_minus_g: array of the u-g color band
    g_mag: array of absolute magnitudes in the g-band
    feh: array of metalicities
    
    
    Returns:
    --------
    kgiant_mask: boolean mask of K-Giant stars
    
    """
    
    kgiant_mask = ((g_minus_r > 0.5) &
                   (g_minus_r < 1.3) &
                   (u_minus_g > 0.5) &
                   (u_minus_g < 3.5) &
                   (g_minus_r <= ((0.087 * feh * 2)
                                  + (0.39 * feh) + 0.96)) & 
                   (g_mag < 4))
    
    return kgiant_mask


