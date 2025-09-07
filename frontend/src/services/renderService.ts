// Interface for render request
export interface RenderRequest {
  input_type: 'SMILES' | 'SMARTS';
  data: string;
}

// Interface for render response
export interface RenderResponse {
  input_type: string;
  data: string;
  image_base64: string;
  success: boolean;
}

class RenderService {
  private baseUrl: string;

  constructor() {
    this.baseUrl = process.env.REACT_APP_BACKEND_URL || 'http://localhost:8000';
  }

  /**
   * Render chemical structures (SMILES or SMARTS) as base64-encoded PNG images
   * 
   * @param request - The render request containing input_type and data
   * @returns Promise<RenderResponse> - The rendered image data
   * 
   * @example
   * // Render SMILES
   * const result = await renderService.renderChemicalStructure({
   *   input_type: 'SMILES',
   *   data: 'CCO'
   * });
   * 
   * @example
   * // Render SMARTS
   * const result = await renderService.renderChemicalStructure({
   *   input_type: 'SMARTS', 
   *   data: '[C:1]=[O:2]>>[C:1][OH:2]'
   * });
   */
  async renderChemicalStructure(request: RenderRequest): Promise<RenderResponse> {
    try {
      // Validate input type
      if (!['SMILES', 'SMARTS'].includes(request.input_type)) {
        throw new Error(`Invalid input_type '${request.input_type}'. Must be 'SMILES' or 'SMARTS'.`);
      }

      // Validate data
      if (!request.data || request.data.trim() === '') {
        throw new Error('Chemical structure data cannot be empty');
      }

      const response = await fetch(`${this.baseUrl}/render`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(request),
      });

      if (!response.ok) {
        const errorText = await response.text();
        let errorMessage = `HTTP error! status: ${response.status}`;
        
        try {
          const errorData = JSON.parse(errorText);
          errorMessage = errorData.detail || errorMessage;
        } catch {
          // If error text is not JSON, use the raw text
          errorMessage = errorText || errorMessage;
        }
        
        throw new Error(errorMessage);
      }

      const result: RenderResponse = await response.json();
      
      // Validate response
      if (!result.success || !result.image_base64) {
        throw new Error('Rendering failed: Invalid response from server');
      }

      return result;
    } catch (error) {
      console.error('Error rendering chemical structure:', error);
      throw error;
    }
  }

  /**
   * Render multiple chemical structures in batch
   * 
   * @param requests - Array of render requests
   * @returns Promise<RenderResponse[]> - Array of rendered results
   */
  async renderMultipleStructures(requests: RenderRequest[]): Promise<RenderResponse[]> {
    const results: RenderResponse[] = [];
    const errors: { index: number; error: string }[] = [];

    for (let i = 0; i < requests.length; i++) {
      try {
        const result = await this.renderChemicalStructure(requests[i]);
        results.push(result);
      } catch (error) {
        const errorMessage = error instanceof Error ? error.message : 'Unknown error';
        errors.push({ index: i, error: errorMessage });
        // Push a failed result to maintain array length
        results.push({
          input_type: requests[i].input_type,
          data: requests[i].data,
          image_base64: '',
          success: false
        });
      }
    }

    if (errors.length > 0) {
      console.warn('Some structures failed to render:', errors);
    }

    return results;
  }

  /**
   * Check if a string is a valid SMILES notation (basic validation)
   * 
   * @param smiles - The SMILES string to validate
   * @returns boolean - Whether the string appears to be valid SMILES
   */
  isValidSMILES(smiles: string): boolean {
    if (!smiles || smiles.trim() === '') return false;
    
    // Basic SMILES validation - contains only valid characters
    const smilesPattern = /^[A-Za-z0-9@+\-\[\]()=#$:.\/\\%]+$/;
    return smilesPattern.test(smiles.trim());
  }

  /**
   * Check if a string is a valid SMARTS notation (basic validation)
   * 
   * @param smarts - The SMARTS string to validate
   * @returns boolean - Whether the string appears to be valid SMARTS
   */
  isValidSMARTS(smarts: string): boolean {
    if (!smarts || smarts.trim() === '') return false;
    
    // Basic SMARTS validation - should contain >> for reactions
    const smartsPattern = /^[A-Za-z0-9@+\-\[\]()=#$:.\/\\%>!&,;]+$/;
    return smartsPattern.test(smarts.trim());
  }

  /**
   * Create a data URL from base64 image data for display
   * 
   * @param base64Data - The base64 encoded image data
   * @param mimeType - The MIME type (default: 'image/png')
   * @returns string - Data URL that can be used in img src
   */
  createImageDataUrl(base64Data: string, mimeType: string = 'image/png'): string {
    // Remove data URL prefix if it exists
    const cleanBase64 = base64Data.replace(/^data:image\/[^;]+;base64,/, '');
    return `data:${mimeType};base64,${cleanBase64}`;
  }
}

// Export a singleton instance
export const renderService = new RenderService();
export default renderService;
