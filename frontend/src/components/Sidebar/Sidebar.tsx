import React, { useState, useEffect } from 'react';
import './Sidebar.css';

interface SidebarProps {
  onMenuSelect?: (menu: string) => void;
  onToggle?: (collapsed: boolean) => void;
  activeMenu?: string;
}

const Sidebar: React.FC<SidebarProps> = ({ onMenuSelect, onToggle, activeMenu: propActiveMenu }) => {
  const [activeMenu, setActiveMenu] = useState<string>(propActiveMenu || 'projects');
  const [isCollapsed, setIsCollapsed] = useState<boolean>(false);

  // Update active menu when prop changes
  useEffect(() => {
    if (propActiveMenu) {
      setActiveMenu(propActiveMenu);
    }
  }, [propActiveMenu]);

  const menuItems = [
    {
      id: 'projects',
      label: 'Projects',
      icon: '▦'
    },
    {
      id: 'dashboard',
      label: 'Dashboard',
      icon: '▢'
    },
    {
      id: 'settings',
      label: 'Settings',
      icon: '⚙'
    }
  ];

  const handleMenuClick = (menuId: string) => {
    setActiveMenu(menuId);
    if (onMenuSelect) {
      onMenuSelect(menuId);
    }
  };

  const toggleSidebar = () => {
    const newCollapsedState = !isCollapsed;
    setIsCollapsed(newCollapsedState);
    if (onToggle) {
      onToggle(newCollapsedState);
    }
  };

  return (
    <div className={`sidebar ${isCollapsed ? 'collapsed' : ''}`}>
      <div className="sidebar-header">
        <h2 className="sidebar-title">{!isCollapsed && 'Reactor'}</h2>
        <button className="sidebar-toggle" onClick={toggleSidebar}>
          {isCollapsed ? '→' : '←'}
        </button>
      </div>
      
      <nav className="sidebar-nav">
        <ul className="sidebar-menu">
          {menuItems.map((item) => (
            <li key={item.id} className="sidebar-menu-item">
              <button
                className={`sidebar-menu-button ${
                  activeMenu === item.id ? 'active' : ''
                }`}
                onClick={() => handleMenuClick(item.id)}
              >
                <span className="sidebar-menu-icon">{item.icon}</span>
                {!isCollapsed && <span className="sidebar-menu-label">{item.label}</span>}
                {isCollapsed && <span className="sidebar-menu-badge">{item.label}</span>}
              </button>
            </li>
          ))}
        </ul>
      </nav>
    </div>
  );
};

export default Sidebar;
