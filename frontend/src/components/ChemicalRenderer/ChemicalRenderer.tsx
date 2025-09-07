import React, { useState, useEffect } from 'react';
import { renderService, type RenderRequest, type RenderResponse } from '../../services';
import './ChemicalRenderer.css';

interface ChemicalRendererProps {
  /** The chemical structure data to render */
  data: string;
  /** Type of input - SMILES or SMARTS */
  inputType: 'SMILES' | 'SMARTS';
  /** Optional width for the rendered image */
  width?: number;
  /** Optional height for the rendered image */
  height?: number;
  /** Optional alt text for the image */
  alt?: string;
  /** Optional CSS class name */
  className?: string;
  /** Callback when rendering completes successfully */
  onRenderSuccess?: (result: RenderResponse) => void;
  /** Callback when rendering fails */
  onRenderError?: (error: string) => void;
}

const ChemicalRenderer: React.FC<ChemicalRendererProps> = ({
  data,
  inputType,
  width,
  height,
  alt,
  className = '',
  onRenderSuccess,
  onRenderError
}) => {
  const [imageUrl, setImageUrl] = useState<string | null>(null);
  const [loading, setLoading] = useState<boolean>(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (!data || data.trim() === '') {
      setError('No chemical structure data provided');
      setLoading(false);
      return;
    }

    const renderStructure = async () => {
      try {
        setLoading(true);
        setError(null);
        setImageUrl(null);

        const request: RenderRequest = {
          input_type: inputType,
          data: data.trim()
        };

        const result = await renderService.renderChemicalStructure(request);
        
        if (result.success && result.image_base64) {
          const dataUrl = renderService.createImageDataUrl(result.image_base64);
          setImageUrl(dataUrl);
          
          if (onRenderSuccess) {
            onRenderSuccess(result);
          }
        } else {
          throw new Error('Rendering failed: No image data received');
        }
      } catch (err) {
        const errorMessage = err instanceof Error ? err.message : 'Unknown rendering error';
        setError(errorMessage);
        
        if (onRenderError) {
          onRenderError(errorMessage);
        }
      } finally {
        setLoading(false);
      }
    };

    renderStructure();
  }, [data, inputType, onRenderSuccess, onRenderError]);

  const containerClasses = `chemical-renderer ${className}`.trim();

  if (loading) {
    return (
      <div className={containerClasses}>
        <div className="chemical-renderer-loading">
          <div className="loading-spinner"></div>
          <p>Rendering {inputType} structure...</p>
        </div>
      </div>
    );
  }

  if (error) {
    return (
      <div className={containerClasses}>
        <div className="chemical-renderer-error">
          <div className="error-icon">⚠️</div>
          <div className="error-content">
            <h4>Rendering Failed</h4>
            <p>{error}</p>
            <div className="error-data">
              <strong>{inputType}:</strong> <code>{data}</code>
            </div>
          </div>
        </div>
      </div>
    );
  }

  if (!imageUrl) {
    return (
      <div className={containerClasses}>
        <div className="chemical-renderer-empty">
          <p>No structure to display</p>
        </div>
      </div>
    );
  }

  return (
    <div className={containerClasses}>
      <div className="chemical-renderer-success">
        <div className="structure-info">
          <span className="structure-type">{inputType}</span>
          <span className="structure-data">{data}</span>
        </div>
        <div className="structure-image">
          <img
            src={imageUrl}
            alt={alt || `${inputType} structure: ${data}`}
            width={width}
            height={height}
            onError={() => setError('Failed to display rendered image')}
          />
        </div>
      </div>
    </div>
  );
};

export default ChemicalRenderer;
